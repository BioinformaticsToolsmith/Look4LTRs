/*
 * Test.cpp
 *
 *  Created on: Nov 29, 2020
 *      Author: Hani Z. Girgis, PhD
 *     Purpose: Covert a multi-sequence alignment produced by clustalo to a Phylip distance matrix
 */

#include <string>
#include <iostream>

#include "../FastaReader.h"
#include "../Matrix.h"

double calculateId(std::string s1, std::string s2, int l) {
	double m = 0.0;
	int a = l;
	for (int i = 0; i < l; i++) {
		if (s1[i] == s2[i]) {
			if (s1[i] == '-') {
				a -= 1;
			} else {
				m += 1;
			}
		}
	}
	return m / a;
}

int main(int argc, char *argv[]) {
	// argv[1]: inputFile
	// argv[2]: outputFile
	if (argc != 3) {
		std::cout << "Use: " << argv[0]
				<< " input_file.fasta output_file.phylip";
		std::cout << std::endl;
		exit(1);
	}

	// Read sequences
	FastaReader reader(std::string(argv[1]), 10000);
	Block *block = reader.read();
	int bSize = block->size();

	// Fill the matrix
	Matrix n(bSize, bSize, 0.0);
	for (int i = 0; i < bSize; i++) {
		if (i % 100 == 0) {
			std::cout << i << std::endl;
		}
		string s1 = *block->at(i).second;
		int l = s1.length();
		for (int j = i + 1; j < bSize; j++) {
			double h = 1.00 - calculateId(s1, *block->at(j).second, l);
			n(i, j) = h;
			n(j, i) = h;
		}
	}

	// Write results
	std::ofstream out(argv[2], std::ios::out);
	out << bSize << std::endl;
	for (int i = 0; i < bSize; i++) {
		out << block->at(i).first->substr(1) << "\t";
		for (int j = 0; j < bSize; j++) {
			out << std::setprecision(8) << fixed << n(i, j);
			if (j != bSize - 1) {
				out << " ";
			}
		}
		out << std::endl;
	}
	out.close();
}
