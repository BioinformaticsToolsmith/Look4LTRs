/*
 * ClustEval.cpp
 *
 *  Created on: Jun 9, 2021
 *      Author: Hani Z. Girgis, PhD
 * This file includes the header and the source!!!
 */

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <thread>

#include "../FastaReader.h"
#include "../meshclust/ClusterInfo.h"
#include "../meshclust/ClusterEvaluator.h"
#include "../meshclust/ClusteringUtil.h"

#include "../GlobAlignE.h"

class ClustEvalAlignment {
private:
	std::string dbFile;
	int threadNum;
	int dataSize;
	Block *data;
	unordered_map<string, int> *table;

	std::vector<ClusterInfo*> makeClusterList(std::string);
	unsigned int countSequenes(std::string);
	double* scoreSeqVsCluster(std::string, ClusterInfo*);

public:
	ClustEvalAlignment(std::string, int);
	virtual ~ClustEvalAlignment();
	void processClusterFile(std::string, std::string);

};

ClustEvalAlignment::ClustEvalAlignment(string database, int cores) {
	dbFile = database;
	threadNum = cores;
	table = new unordered_map<string, int>();

	// Read database
	cout << "Reading database ..." << endl;
	dataSize = countSequenes(dbFile);
	FastaReader dbReader(dbFile, dataSize);
	data = dbReader.read();
	int blockSize = data->size();
	if (blockSize != dataSize) {
		std::cerr << "ClustEvalAlignment: Something wrong while reading data.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// Build the hash table: header -> index
	cout << "Building table ..." << endl;
	for (int i = 0; i < blockSize; i++) {
		table->emplace(*data->at(i).first, i);
	}
}

ClustEvalAlignment::~ClustEvalAlignment() {
	// do nothing
}

unsigned int ClustEvalAlignment::countSequenes(std::string file) {
	std::ifstream in(file.c_str());
	unsigned int count = 0;
	std::string line;
	while (getline(in, line)) {
		if (line[0] == '>') {
			++count;
		}
	}
	return count;
}

std::vector<ClusterInfo*> ClustEvalAlignment::makeClusterList(
		std::string fileName) {
	std::vector<ClusterInfo*> list;
	ifstream in(fileName);
	int identifier;
	string info;
	double score;
	char membership;
	while (in >> identifier >> info >> score >> membership) {
		if (identifier > list.size()) {
			list.push_back(new ClusterInfo(identifier));
		}
		auto cluster = list.at(identifier - 1);
		cluster->addMember(new string(info), membership == 'C' ? true : false,
				0.0, ClusteringUtil::MEMBER);
	}
	in.close();

	return list;
}

/**
 * The client has to free the memory of the returned array
 */
double* ClustEvalAlignment::scoreSeqVsCluster(std::string header,
		ClusterInfo *cluster) {

	// Find the index of the header of interest
	auto itr = table->find(header);
	if (itr == table->end()) {
		cerr << "Center " << header << " is not in the table." << endl;
		throw std::exception();
	}
	string *seq1 = data->at(itr->second).second;

	// Collect histograms of members
	int memberNum = cluster->getSize();
	auto memberList = cluster->getMemberList();
	double *scoreList = new double[memberNum];
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < memberNum; i++) {
		std::string member = *memberList[i];
		auto itr = table->find(member);
		if (itr == table->end()) {
			cerr << "Sequence " << member << " is not in the table." << endl;
			throw std::exception();
		}
		string *seq2 = data->at(itr->second).second;

		GlobAlignE align(seq1->c_str(), 0, seq1->size() - 1, seq2->c_str(), 0,
				seq2->size() - 1, 1, -1, 4, 1);
		scoreList[i] = align.getIdentity();
	}

	return scoreList;
}

void ClustEvalAlignment::processClusterFile(std::string clusterFile,
		std::string outFile) {
	// Align all centers versus all centers
	std::vector<ClusterInfo*> clusterList = makeClusterList(clusterFile);
	std::vector<std::string*> centerList(clusterList.size(), nullptr);
	int clusterNum = clusterList.size();

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < clusterNum; i++) {
		std::string header = *clusterList[i]->getCenter();
		auto itr = table->find(header);
		if (itr == table->end()) {
			cerr << "Center " << header << " is not in the table." << endl;
			throw std::exception();
		}
		centerList[i] = data->at(itr->second).second;
	}

	// Initialize the all vs. all matrix
	Matrix cvc(clusterNum, clusterNum, 0.0);
	for (int i = 0; i < clusterNum; i++) {
		cvc(i, i) = 1.0;
	}

	// Score all vs. all
	for (int i = 0; i < clusterNum; i++) {
		std::string *seq1 = centerList[i];
		int s1 = 0;
		int e1 = seq1->size() - 1;
#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int j = i + 1; j < clusterNum; j++) {
			std::string *seq2 = centerList[j];
			GlobAlignE align(seq1->c_str(), s1, e1, seq2->c_str(), 0,
					seq2->size() - 1, 1, -1, 4, 1);
			cvc(i, j) = cvc(j, i) = align.getIdentity();
		}
	}

	// Score members versus their center
	for (auto clusterPtr : clusterList) {
		double *scoreList = scoreSeqVsCluster(*clusterPtr->getCenter(),
				clusterPtr);
		clusterPtr->updateScoreWithCenter(scoreList, clusterPtr->getSize());
		delete[] scoreList;
	}

	// Score members versus closest neighbor
	for (auto clusterPtr : clusterList) {
		int i = clusterPtr->getIdentifier() - 1;
		double max = -1.0;
		int index = -1;
		for (int j = 0; j < clusterNum; j++) {
			if (i != j) {
				if (cvc(i, j) > max) {
					max = cvc(i, j);
					index = j;
				}
			}
		}

		if (max == -1.0 || index == -1) {
			std::cerr << "Could not find the closest neighbor." << std::endl;
			throw std::exception();
		}

		string neighbor = *clusterList.at(index)->getCenter();
		double *scoreList = scoreSeqVsCluster(neighbor, clusterPtr);
		clusterPtr->updateScoreWithNeighbor(scoreList, clusterPtr->getSize());
		delete[] scoreList;
	}

	// Evaluate
	ClusterEvaluator evaluator(cvc, clusterList, dataSize);
	double db = evaluator.daviesBouldin();
	double dn = evaluator.dunn();
	double sil = evaluator.silhoutte();
	double ratio = evaluator.clusterRatio();
	double intraValue = evaluator.intra();
	double interValue = evaluator.inter();

	// Write results
	ofstream out(outFile);
	out << std::setprecision(4);
	out << "Davies-Bouldin\t" << db << std::endl;
	out << "Dunn\t" << dn << std::endl;
	out << "Silhoutte\t" << sil << std::endl;
	out << "Coverage\t" << ratio << std::endl;
	out << "Intra\t" << intraValue << std::endl;
	out << "Inter\t" << interValue << std::endl;
	out.close();
}

int main(int argc, char *argv[]) {
	string dbFile("");
	string clusterFile("");
	string modelFile("");
	string outFile("");
	int cores = std::thread::hardware_concurrency();

	if (argc != 9) {
		cerr << "Usage: " << argv[0];
		cerr << " -d database_file -g cluster_file -o output_file -c core_number";
		cerr << endl;
		exit(1);
	}

	for (int i = 1; i < argc; i += 2) {
		switch (argv[i][1]) {
		case 'd': {
			dbFile = std::string(argv[i + 1]);
		}
			break;
		case 'g': {
			clusterFile = std::string(argv[i + 1]);
		}
			break;

		case 'o': {
			outFile = std::string(argv[i + 1]);
		}
			break;

		case 'c': {
			cores = atoi(argv[i + 1]);
		}
			break;
		}
	}

	ClustEvalAlignment evaluator(dbFile, cores);
	evaluator.processClusterFile(clusterFile, outFile);

	return 0;
}
