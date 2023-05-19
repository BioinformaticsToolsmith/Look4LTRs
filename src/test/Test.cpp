/*
 * Test.cpp
 *
 *  Created on: Sep 5, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#include <string>
#include <iostream>

#include "../FastaReader.h"
#include "../KmerHistogram.h"
#include "../Matrix.h"

int main(int argc, char *argv[]) {
	int listSize = 46340; // Works up to 46340;
	int threadNum = 16;

	std::cout << "About allocating the matrix ..." << std::endl;
	Matrix m(listSize, listSize, 1.0);
	std::cout << "Done allocating the matrix ..." << std::endl;

	m = Matrix();

	while (true)
		;

//	for (int i = 0; i < listSize; i++) {
//
////#pragma omp parallel for schedule(static) num_threads(threadNum)
//		for (int j = i + 1; j < listSize; j++) {
//			m->at(i, j) = 0.0;
//			m->at(j, i) = 0.0;
//		}
//	}

//	FastaReader reader(
//			std::string(
//					"/home/zakaria/Data/Identity/Keratin/keratin_small.fasta"),
//			1);
//	Block *block = reader.read();
//	int k = 2;
//	KmerHistogram<int, int> table(k);
//	const std::string *seqPtr = block->at(0).second;
//	int lastIndex = seqPtr->size() - k + 1;
//	vector<int> hashList;
//	table.hash(seqPtr, 0, lastIndex, &hashList);
//
//	for (int i = 0; i < hashList.size() - 1; i++) {
//		auto item = hashList.at(i);
//		if (item == 2) {
//			std::cout << item << " " << hashList.at(i + 1) << std::endl;
//		}
//	}
}
