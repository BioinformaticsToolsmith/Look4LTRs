/*
 MeShClust v3.0 clusters sequences using the mean shift algorithm and alignment-free identity scores.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * GMM.cpp
 *
 *  Created on: Oct 15, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#include "GMM.h"

GMM::GMM(const vector<double> &d, int n) :
		dataList(d), k(n) {

	if (k < 1) {
		cerr << "The number of clusters cannot be less than 1.";
		cerr << endl;
		throw std::exception();
	}

	double min = *min_element(d.begin(), d.end());
	double max = *max_element(d.begin(), d.end());
	// Place initial centers uniformly in this interval
	double h = (max - min) / (k + 1);

	for (int i = 0; i < k; i++) {
		muList.push_back(min + ((i + 1) * h));
		sigmaList.push_back(h);
	}
	phiList = vector<double>(muList.size(), 1.0 / k);
	gammaMat = Matrix(dataList.size(), k, 1.0 / k);

	run();
}

GMM::GMM(const vector<double> &d, vector<double> &m, vector<double> &s,
		vector<double> &f) :
		dataList(d) {

	if (m.size() != s.size()) {
		cerr << "The number of means does not equal the number of std's.";
		cerr << endl;
		throw std::exception();
	}

	k = m.size();
	muList = m;
	sigmaList = s;
	phiList = f;
	gammaMat = Matrix(dataList.size(), k, 1.0 / k);

	run();
}

GMM::~GMM() {
	// TODO Auto-generated destructor stub
}

void GMM::run() {
	int i = 0;
	for (; i < itr; i++) {
		oldMuList = muList;
		updatePartitions();
		updateParameters();
		if (isConverged()) {
			break;
		}
	}
	cout << "Converged after " << i << " iterations." << endl;
}

/**
 * Update responsibilities
 */
void GMM::updatePartitions() {
	int size = dataList.size();

#pragma omp parallel for schedule(static)
	for (int i = 0; i < size; i++) {
		double sum = 0.0;
		for (int j = 0; j < k; j++) {
			double pj = phiList[j]
					* calcP(dataList[i], muList[j], sigmaList[j]);
			sum += pj;
			gammaMat(i, j) = pj;
		}

		for (int j = 0; j < k; j++) {
			gammaMat(i, j) /= sum;
		}
	}
}

/**
 * Update means, standard deviations, and cluster weights
 */
void GMM::updateParameters() {
	int size = dataList.size();

#pragma omp parallel for schedule(static)
	for (int j = 0; j < k; j++) {
		double gammaSum = 0.0;
		for (int i = 0; i < size; i++) {
			gammaSum += gammaMat(i, j);
		}

		double mu = 0.0;
		double sigma = 0.0;
		for (int i = 0; i < size; i++) {
			mu += gammaMat(i, j) * dataList[i];
			double t = dataList[i] - muList[j];
			sigma += gammaMat(i, j) * t * t;
		}
		muList[j] = mu / gammaSum;
		sigmaList[j] = sqrt(sigma / gammaSum);

		phiList[j] = gammaSum / size;

		if (phiList[j] < 0 || phiList[j] > 1) {
			cerr << "Incorrect cluster weight of " << phiList[j];
			cerr << endl;
			throw std::exception();
		}
	}
}

/**
 * Compare the new means to the old ones.
 */
bool GMM::isConverged() {
	bool r = true;
	for (int i = 0; i < k; i++) {
		if (!Util::isEqual(muList[i], oldMuList[i])) {
			r = false;
			break;
		}
	}
	return r;
}

/**
 * BIC = log(m) x parameter_num - 2 log(likelihood)
 */
double GMM::calcBIC() {
	int m = dataList.size();

	vector<double> list(m, 0.0);
#pragma omp parallel for schedule(static)
	for (int i = 0; i < m; i++) {
		double p = 0.0;
		for (int j = 0; j < k; j++) {
			p += phiList[j] * calcP(dataList[i], muList[j], sigmaList[j]);
		}
		list[i] = log(p);
	}

	double l = 0.0;
	for (int i = 0; i < m; i++) {
		l += list[i];
	}
	// Number of parameters: k means, k standard deviations, k-1 weights
	double n = k + k + k - 1;
	return log(m) * n - 2 * l;
}

vector<double> GMM::getMuList() {
	return muList;
}

vector<double> GMM::getSigmaList() {
	return sigmaList;
}

vector<double> GMM::getPhiList() {
	return phiList;
}

Matrix GMM::getGammaMat() {
	return gammaMat;
}

void GMM::print() {
	for (int i = 0; i < k; i++) {
		cout << "Cluster " << i << ": " << phiList[i] << " " << muList[i] << " "
				<< sigmaList[i] << endl;
	}
	cout << "BIC: " << calcBIC() << endl;
}
