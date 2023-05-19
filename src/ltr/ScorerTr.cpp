/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * ScorerTr.cpp
 *
 *  Created on: Oct 4, 2022
 *      Author: Anthony B. Garza.
 *      Author: Hani Z. Girgis
 *
 * Academic use: Affero General Public License version 1.
 *
 * Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 *
 * Copyright (C) 2022 by the authors.
 */
#include "ScorerTr.h"

ScorerTr::ScorerTr(std::string &seqIn, int motifSizeIn, int minIn, int maxIn) : seq(seqIn)
{
	// Pre-conditions
	assert(motifSizeIn >= 1);
	assert(seqIn.size() >= motifSizeIn);
	assert(minIn <= maxIn);
	assert(minIn >= 1);
	assert(maxIn >= 1);

	// seq = seqIn;
	k = motifSizeIn;
	min = minIn;
	max = maxIn;

	kmerTable = new KmerHistogram<int, int>(k);

	// std::cout << LtrParameters::INITIAL_VALUE << std::endl;
	// std::cout << LtrParameters::INITIAL_SCORE << std::endl;

	int init_score = LtrParameters::INITIAL_SCORE;

	forwardList = new std::vector<int>(seq.size(), init_score);
	backwardList = new std::vector<int>(seq.size(), init_score);

	score();
}

ScorerTr::~ScorerTr()
{
	delete kmerTable;

	forwardList->clear();
	delete forwardList;

	backwardList->clear();
	delete backwardList;
}

void ScorerTr::score()
{
	int *indexList = kmerTable->initialize(LtrParameters::INITIAL_VALUE);
	vector<pair<int, int>> segmentList = kmerTable->makeSegments(&seq);

	for (auto segment : segmentList)
	{
		int start = segment.first;
		int end = segment.second + 1; // End is inclusive, add 1

		if (end - start < k)
		{ // If start = 0, end = 2, and k = 3... the segment is smaller than k and thus should be skipped
			continue;
		}

		vector<int> *hashList = new vector<int>();

		hashList->reserve(end - start - (k - 1));

		//kmerTable->hash(&seq, start, end - k + 2, hashList);
		kmerTable->hash(&seq, start, end - k, hashList);

		for (int i = start; i <= end - k; i++)
		{ // Any position at end - k + 2 or greater can't have a kmer of size k
			int keyHash = hashList->at(i - start);
			// Where is the closest previous copy
			int lastIndex = indexList[keyHash];

			// There is a copy of this kmer and its location is at lastIndex
			if (lastIndex != LtrParameters::INITIAL_VALUE)
			{
				int dist = i - lastIndex;
				if (dist >= min && dist <= max)
				{
					// Look forward
					forwardList->at(lastIndex) = dist;

					// Look backward
					backwardList->at(i) = dist;
				}
			}

			indexList[keyHash] = i;
		}

		hashList->clear();
		delete hashList;
	}
}


vector<int> *ScorerTr::getForwardScores()
{
	return forwardList;
}

vector<int> *ScorerTr::getBackwardScores()
{
	return backwardList;
}

void ScorerTr::printForwardScores(std::string fileName)
{
	std::stringstream s(std::ios_base::out);
	for (auto score : *forwardList) {
		s << score << std::endl;
	}


	ofstream out(fileName);
	out.write(s.str().c_str(), s.str().length());
	out.close();
}

void ScorerTr::printBackwardScores(std::string fileName)
{
	std::stringstream s(std::ios_base::out);
	for (auto score : *backwardList) {
		s << score << std::endl;
	}


	ofstream out(fileName);
	out.write(s.str().c_str(), s.str().length());
	out.close();
}
