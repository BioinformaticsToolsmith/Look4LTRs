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
 * GMM.h
 *
 *  Created on: Oct 15, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef SRC_MESHCLUST_GMM_H_
#define SRC_MESHCLUST_GMM_H_

#include <vector>
#include <math.h>
#include <algorithm>

#include "../Matrix.h"
#include "../Util.h"

using namespace std;

class GMM {
private:
	const vector<double> &dataList;
	vector<double> muList;
	vector<double> oldMuList;
	vector<double> sigmaList;
	vector<double> phiList;
	Matrix gammaMat;
	int k;

	const int itr = 200;
	const double sqrtTwoPi = sqrt(2 * M_PI);

	inline double calcP(double x, double mu, double sigma) {
		double xMinusMuSquared = x - mu;
		xMinusMuSquared *= xMinusMuSquared;
		return (1.0 / (sigma * sqrtTwoPi))
				* exp(-xMinusMuSquared / (2 * sigma * sigma));
	}

	void updatePartitions();
	void updateParameters();
	bool isConverged();
	void run();

public:
	GMM(const vector<double>&, int);
	GMM(const vector<double>&, vector<double>&, vector<double>&,
			vector<double>&);
	virtual ~GMM();
	double calcBIC();
	vector<double> getMuList();
	vector<double> getSigmaList();
	vector<double> getPhiList();
	Matrix getGammaMat();
	void print();
};

#endif /* SRC_MESHCLUST_GMM_H_ */
