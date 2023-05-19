/*
 * ClustEval.cpp
 *
 *  Created on: Jun 9, 2021
 *      Author: Hani Z. Girgis, PhD
 * This file includes the header and the source!!!
 */

// ToDo:
// Remove small cluster < 5 members or singles
//
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <thread>

#include "../FastaReader.h"
#include "../IdentityCalculator.h"
#include "../Serializer.h"
#include "../meshclust/ClusterInfo.h"
#include "../meshclust/ClusterEvaluator.h"
#include "../meshclust/ClusteringUtil.h"

template<class V>
class ClustEval {
private:
	V **kHistList;
	uint64_t **monoHistList;
	std::string **infoList;
	int *lenList;

	std::string dbFile;
	int threadNum;

	int dataSize;

	IdentityCalculator<V> *identity;
	unordered_map<string, int> *table;

	std::vector<ClusterInfo*> makeClusterList(std::string);
	unsigned int countSequenes(std::string);
	double* scoreSeqVsCluster(std::string, ClusterInfo*);

public:
	ClustEval(Serializer&, std::string, int);
	virtual ~ClustEval();
	void processClusterFile(std::string, std::string);

};

template<class V>
ClustEval<V>::ClustEval(Serializer &s, string database, int cores) {
	dbFile = database;
	threadNum = cores;

	// No threshold is provided, cannot skip, and cannot relax
	// Why no threshold is provided? Because we need not to skip.
	identity = new IdentityCalculator<V>(s, 0.0, false, false);

	table = new unordered_map<string, int>();

	// Read database
	cout << "Reading database ..." << endl;
	dataSize = countSequenes(dbFile);
	FastaReader dbReader(dbFile, dataSize);
	Block *block = dbReader.read();
	int blockSize = block->size();
	if (blockSize != dataSize) {
		std::cerr << "ClustEval: Something wrong while reading data.";
		std::cerr << std::endl;
		throw std::exception();
	}

	auto tuple = identity->unpackBlock(block, cores);
	kHistList = get<0>(tuple);
	monoHistList = get<1>(tuple);
	infoList = get<2>(tuple);
	lenList = get<3>(tuple);

	// Build the hash table: header -> index
	cout << "Building table ..." << endl;
	for (int i = 0; i < blockSize; i++) {
		table->emplace(*infoList[i], i);
	}
}

template<class V>
ClustEval<V>::~ClustEval() {
	delete identity;
}

template<class V>
unsigned int ClustEval<V>::countSequenes(std::string file) {
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

template<class V>
std::vector<ClusterInfo*> ClustEval<V>::makeClusterList(std::string fileName) {
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
template<class V>
double* ClustEval<V>::scoreSeqVsCluster(std::string header,
		ClusterInfo *cluster) {

	// Find the index of the header of interest
	auto itr = table->find(header);
	if (itr == table->end()) {
		cerr << "Center " << header << " is not in the table." << endl;
		throw std::exception();
	}
	int c = itr->second;

	// Collect histograms of members
	int memberNum = cluster->getSize();
	V *kHistListCluster[memberNum];
	uint64_t *monoHistListCluster[memberNum];
	std::string *infoListCluster[memberNum];
	int lenListCluster[memberNum];

	auto memberList = cluster->getMemberList();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < memberNum; i++) {
		std::string member = *memberList[i];
		auto itr = table->find(member);
		if (itr == table->end()) {
			cerr << "Sequence " << member << " is not in the table." << endl;
			throw std::exception();
		}
		int j = itr->second;

		kHistListCluster[i] = kHistList[j];
		monoHistListCluster[i] = monoHistList[j];
		infoListCluster[i] = infoList[j];
		lenListCluster[i] = lenList[j];
	}

	double *scoreList = identity->score(kHistList[c], kHistListCluster,
			monoHistList[c], monoHistListCluster, memberNum, threadNum,
			lenList[c], lenListCluster);

	return scoreList;
}

template<class V>
void ClustEval<V>::processClusterFile(std::string clusterFile,
		std::string outFile) {
	// Align all centers versus all centers
	std::vector<ClusterInfo*> clusterList = makeClusterList(clusterFile);
	int clusterNum = clusterList.size();
	V *kHistListCVC[clusterNum];
	uint64_t *monoHistListCVC[clusterNum];
	std::string *infoListCVC[clusterNum];
	int lenListCVC[clusterNum];
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < clusterNum; i++) {
		std::string header = *clusterList[i]->getCenter();
		auto itr = table->find(header);
		if (itr == table->end()) {
			cerr << "Center " << header << " is not in the table." << endl;
			throw std::exception();
		}
		int j = itr->second;
		kHistListCVC[i] = kHistList[j];
		monoHistListCVC[i] = monoHistList[j];
		infoListCVC[i] = infoList[j];
		lenListCVC[i] = lenList[j];
	}
	Matrix cvc = identity->score(kHistListCVC, monoHistListCVC, clusterNum,
			threadNum, lenListCVC);

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

	if (argc != 11) {
		cerr << "Usage: " << argv[0];
		cerr
				<< " -d database_file -g cluster_file -m model_file -o output_file -c core_number";
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

		case 'm': {
			modelFile = std::string(argv[i + 1]);
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

	// Build model
	Serializer s(modelFile);
	int64_t maxLength = s.getMaxLength();

	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		ClustEval<int8_t> evaluator(s, dbFile, cores);
		evaluator.processClusterFile(clusterFile, outFile);
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		ClustEval<int16_t> evaluator(s, dbFile, cores);
		evaluator.processClusterFile(clusterFile, outFile);
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		ClustEval<int32_t> evaluator(s, dbFile, cores);
		evaluator.processClusterFile(clusterFile, outFile);
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		ClustEval<int64_t> evaluator(s, dbFile, cores);
		evaluator.processClusterFile(clusterFile, outFile);
	} else {
		std::cout << "Warning: Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		ClustEval<int64_t> evaluator(s, dbFile, cores);
		evaluator.processClusterFile(clusterFile, outFile);
	}
	return 0;
}
