/*
 * GenerateClusters.cpp
 *
 *  Created on: May 29, 2021
 *      Author: Dr. Hani Z. Girgis, the Bioinformatics Toolsmith Laboratory, Texas A&M University-Kingsville
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <atomic>
#include <chrono>

#include "../Mutator.h"
#include "../GlobAlignE.h"

using namespace std;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine gen(seed);
uniform_real_distribution<> zeroOneRand;
// A list of the generated templates
vector<string> templateList;
// A list of cluster sizes
vector<int> clustSizeList;
vector<int> templateSizeList;

/**
 * Generate a random nucleotide
 */
char getRandomNucleotide() {
	double p = zeroOneRand(gen);
	char r;
	if (p <= 0.25) {
		r = 'A';
	} else if (p <= 0.50) {
		r = 'C';
	} else if (p <= 0.75) {
		r = 'G';
	} else {
		r = 'T';
	}
	return r;
}

/**
 * Generate a random template
 */
string makeRandomTemplate(int len) {
	string seq("");
	for (int i = 0; i < len; i++) {
		seq.append(1, getRandomNucleotide());
	}
	return seq;
}

/**
 * Make sure that a template is different from the ones we have so far
 */
bool isValidTemplate(string s, double idThreshold) {
	bool r = true;
	std::atomic<bool> canContinue = false;

#pragma omp parallel for schedule(static)
	for (int i = 0; i < templateList.size(); i++) {
		if (canContinue) {
			continue;
		}

		string t = templateList[i];
		GlobAlignE align(s.c_str(), 0, s.size() - 1, t.c_str(), 0, t.size() - 1,
				1, -1, 4, 1);

		if (align.getIdentity() >= idThreshold) {
			// cerr << "Excluding an invalid template: " << align.getIdentity() << " " << idThreshold << endl;
			r = false;
			canContinue = true;
		}

	}
	return r;
}

int mean(vector<int> list) {
	double average = 0.0;
	for (int x : list) {
		average += x;
	}
	return round(average / list.size());
}

int main(int argc, char *argv[]) {
	if (argc < 9) {
		cerr
				<< "The program needs min_sequence_length (>=1) max_sequence_length (>=1) cluster_num (>=1) min_cluster_size (>=1) max_cluster_size (>=1) max_mutation_rate (0-1) include_template (y or n) output_file (fasta)";
		cerr << endl;
		exit(1);
	}

	// minimum template length
	int minSeqLen = atoi(argv[1]);

	// maximum template length
	int maxSeqLen = atoi(argv[2]);

	// Number of clusters
	int clusterNum = atoi(argv[3]);
	// Minimum cluster size
	int minClustSize = atoi(argv[4]);
	// Maximum cluster size
	int maxClustSize = atoi(argv[5]);
	// Maximum mutation rate
	double mutRate = atof(argv[6]);
	// If true, the original template is included in the generated cluster as a member
	bool canIncludeTemplate = argv[7][0] == 'y' ? true : false;
	// Output file including the generated clusters
	string outputFile(argv[8]);

	// Validate arguments
	if (minSeqLen <= 0) {
		cerr << "Minimum sequence length must be 1 or more." << endl;
		exit(1);
	}

	if (maxSeqLen <= 0) {
		cerr << "maximum sequence length must be 1 or more." << endl;
		exit(1);
	}

	if (minSeqLen > maxSeqLen) {
		cerr << "Maximum sequence length cannot be smaller than its minimum.";
		cerr << endl;
		exit(1);
	}

	if (clusterNum <= 0) {
		cerr << "Cluster number must be 1 or more." << endl;
		exit(1);
	}

	if (minClustSize <= 0) {
		cerr << "Minimum cluster size must be 1 or more." << endl;
		exit(1);
	}

	if (maxClustSize <= 0) {
		cerr << "Maximum cluster size must be 1 or more." << endl;
		exit(1);
	}

	if (maxClustSize < minClustSize) {
		cerr << "Maximum cluster size cannot be smaller than its minimum.";
		cerr << endl;
		exit(1);
	}

	if (mutRate < 0.0 || mutRate > 1.0) {
		cerr << "Mutation rate must be between zero and one.";
		cerr << endl;
		exit(1);
	}

	if (!(argv[7][0] == 'y' || argv[7][0] == 'n')) {
		cerr << "To include the template, provide y or n.";
		cerr << endl;
		exit(1);
	}

	// Fill the list of templates
	double idThreshold = 1.0 - mutRate - 0.1; //0.05; // Make sure that the cluster are separated enough

	uniform_int_distribution<int> templateRand(minSeqLen, maxSeqLen);

	for (int i = 0; i < clusterNum; i++) {
		int seqLen = templateRand(gen);
		//cout << "Current template length: " << seqLen << endl;
		while (true) {
			string s = makeRandomTemplate(seqLen);
			if (isValidTemplate(s, idThreshold)) {
				templateList.push_back(s);
				templateSizeList.push_back(seqLen);
				break;
			} else {
				// cerr << "Excluding an invalid template." << endl;
				// cerr << ".";
			}
		}
	}

	// Open output file
	ofstream out(outputFile);

	int total = 0;
	for (int i = 0; i < templateList.size(); i++) {
		string seq = templateList[i];
		if (canIncludeTemplate) {
			out << ">template_" << i + 1 << endl;
			out << seq << endl;
		}

		Mutator mutator(&seq, 5, i, 2);
		mutator.enableSinglePoint();
		mutator.enableBlock();
		uniform_int_distribution<int> minMaxRand(minClustSize, maxClustSize);
		int clustSize = minMaxRand(gen);
		clustSizeList.push_back(clustSize);

		uniform_real_distribution<> mutRateRand(0.009, mutRate);
		for (int j = 0; j < clustSize;) {
			pair<string*, double> seqIdentity = mutator.mutateSequence(
					mutRateRand(gen));
			string *mutSeq = seqIdentity.first;
			double identity = seqIdentity.second;
			if (1.0 - identity <= mutRate) {
				out << ">template_" << i + 1 << "_" << j + 1 /*<< " " << identity*/
				<< endl;
				out << *mutSeq << endl;
				j++;
				total++;
			}
			delete mutSeq;
		}
	}

	// Clean up
	out.close();

	cout << "Average template length: " << mean(templateSizeList) << endl;
	cout << "Minimum template length: "
			<< *std::min_element(templateSizeList.begin(),
					templateSizeList.end()) << endl;
	cout << "Maximum template length: "
			<< *std::max_element(templateSizeList.begin(),
					templateSizeList.end()) << endl;

	cout << "Average cluster size: " << mean(clustSizeList) << endl;
	cout << "Minimum cluster size: "
			<< *std::min_element(clustSizeList.begin(), clustSizeList.end())
			<< endl;
	cout << "Maximum cluster size: "
			<< *std::max_element(clustSizeList.begin(), clustSizeList.end())
			<< endl;
	cout << "Cluster number: " << clusterNum << endl;
	cout << "Total number of sequences: " << total << endl;
	// Some of the information displayed for testing
}
