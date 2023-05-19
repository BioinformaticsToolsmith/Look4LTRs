/*
 * makeClusters.cpp
 *
 *  Created on: Mar 6, 2021
 *      Author: Hani Z. Girgis, PhD
 *     Purpose: Clusters the output of Identity by finding connected components.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include "../Matrix.h"
#include "../meshclust/ClusteringUtil.h"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 4) {
		cerr << "The program needs a threshold score and ";
		cerr << "an all_vs_all_file output_file." << endl;
		cerr << "Example: makeclusters 0.75 allVsAll.txt result.txt" << endl;
		cerr << "The file does not need to include all possible pairs." << endl;
		exit(1);
	}

	double threshold = atof(argv[1]);
	ifstream fin(argv[2]);

	// Count unique sequences that meet the criterion (score >= threshold)
	cout << "Count unique sequences with scores >= threshold ..." << endl;
	string id1, id2;
	double score;
	unordered_map<string, double> label_dict;
	int seq_index = 0;
	while (fin >> id1 >> id2 >> score) {
		//if (score >= threshold) {
		if (label_dict.find(id1) == label_dict.end()) {
			label_dict[id1] = seq_index++;
		}

		if (label_dict.find(id2) == label_dict.end()) {
			label_dict[id2] = seq_index++;
		}
		//}
	}

	// Build the matrix
	int size = label_dict.size();
	cout << "Fill matrix " << size << " x " << size << " ..." << endl;
	Matrix m(size, size);
	fin.clear();
	fin.seekg(0);

	while (fin >> id1 >> id2 >> score) {
		//if (score >= threshold) {
		int index1 = label_dict[id1];
		int index2 = label_dict[id2];
		if (index1 < index2) {
			m.at(index1, index2) = score;
		} else {
			m.at(index2, index1) = score;
		}
		//}
	}
	fin.close();

	// Find connected component
	cout << "Cluster ..." << endl;
	auto p = ClusteringUtil::findConnectedComponents(m, threshold);
	vector<int> flagList = p.first;
	int compNum = p.second;

	// Make a list of labels (sequence names) according to the matrix order
	vector<string> label_list(size, "");
	unordered_map<string, double>::iterator iter;
	for (iter = label_dict.begin(); iter != label_dict.end(); iter++) {
		label_list[iter->second] = iter->first;
	}

	// Write results
	cout << "Write results ..." << endl;
	ofstream fout(argv[3]);
	for (int c = 1; c <= compNum; c++) {
		fout << "Cluster " << c << ":" << endl;
		for (int i = 0; i < size; i++) {
			if (flagList[i] == c) {
				fout << "\t" << label_list[i] << endl;
			}
		}
		fout << endl;
	}

	// Free resources
	cout << "Done." << endl;
	fout.close();
	return 0;
}
