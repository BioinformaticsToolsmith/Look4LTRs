/*
 * IdentityForClustEval.cpp
 *
 *  Created on: Jun 7, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <thread>

#include "../FastaReader.h"
#include "../IdentityCalculator.h"
#include "../Serializer.h"

using namespace std;

unsigned int countSequenes(std::string file) {
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
void align(Serializer &s, string dbFile, string pairFile, string outFile,
		int cores) {
	IdentityCalculator<V> identity(s, 0.0, false, false);

	// Read database
	cout << "Reading database ..." << endl;
	FastaReader dbReader(dbFile, countSequenes(dbFile));
	Block *block = dbReader.read();
	int blockSize = block->size();
	auto tuple = identity.unpackBlock(block, cores);
	V **kHistList = get<0>(tuple);
	uint64_t **monoHistList = get<1>(tuple);
	std::string **infoList = get<2>(tuple);
	int *lenList = get<3>(tuple);

	// Build the hash table: header -> index
	cout << "Building table ..." << endl;
	unordered_map<string, int> table;
	for (int i = 0; i < blockSize; i++) {
		table[*infoList[i]] = i;
	}

	// Read pairs
	// Assumes that a sequence name does not include spaces
	cout << "Reading pairs ..." << endl;
	vector<pair<string, string> > pairList;
	ifstream in(pairFile);
	string first, second;
	while (in >> first >> second) {
		pairList.push_back(make_pair(first, second));
	}
	in.close();

	// Align
	cout << "Calculating identity scores ..." << endl;
	int pairCount = pairList.size();
	vector<double> scoreList(pairCount);
#pragma omp parallel for schedule(static) num_threads(cores)
	for (int i = 0; i < pairCount; i++) {
		auto p = pairList[i];

		auto itr1 = table.find(p.first);
		if (itr1 == table.end()) {
			cerr << "Sequence " << p.first << " is not in the table." << endl;
			throw std::exception();
		}
		auto itr2 = table.find(p.second);
		if (itr2 == table.end()) {
			cerr << "Sequence " << p.second << " is not in the table." << endl;
			throw std::exception();
		}

		int index1 = itr1->second;
		int index2 = itr2->second;

		double ratio = identity.calcRatio(lenList[index1], lenList[index2]);
		scoreList[i] = identity.score(kHistList[index1], kHistList[index2],
				monoHistList[index1], monoHistList[index2], ratio);
	}

	// Write output results
	cout << "Writing results ..." << endl;
	ofstream out(outFile);
	for (int i = 0; i < pairCount; i++) {
		auto p = pairList[i];
		out << p.first << "\t" << p.second << "\t" << scoreList[i] << endl;
	}
	out.close();
}

int main(int argc, char *argv[]) {
	string dbFile("");
	string pairFile("");
	string modelFile("");
	string outFile("");
	int cores = std::thread::hardware_concurrency();

	if (argc != 11) {
		cerr << "Usage: " << argv[0];
		cerr
				<< " -d database_file -p pair_file -m model_file -o output_file -c core_number";
		cerr << endl;
		exit(1);
	}

	for (int i = 1; i < argc; i += 2) {
		switch (argv[i][1]) {
		case 'd': {
			dbFile = std::string(argv[i + 1]);
		}
			break;
		case 'p': {
			pairFile = std::string(argv[i + 1]);
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
		align<int8_t>(s, dbFile, pairFile, outFile, cores);
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		align<int16_t>(s, dbFile, pairFile, outFile, cores);
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		align<int32_t>(s, dbFile, pairFile, outFile, cores);
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		align<int64_t>(s, dbFile, pairFile, outFile, cores);
	} else {
		std::cout << "Warning: Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		align<int64_t>(s, dbFile, pairFile, outFile, cores);
	}
	return 0;
}
