/*
 * ScorerTr.h
 *
 *  Created on: Nov 30, 2012
 *      Author: Hani Zakaria Girgis, PhD
 *      Author: Mr. Anthony B. Garza
 */

#ifndef SCORERTR_H_
#define SCORERTR_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <sstream>

#include "LtrParameters.h"
#include "../KmerHistogram.h"

class ScorerTr
{
private:
	/**
	 * Variables
	 */

	// Reference to sequence to score
	std::string &seq;
	// Size of k-mer
	int k;
	// Minimum score possible; inclusive
	int min;
	// Maximum score possible; inclusive
	int max;

	KmerHistogram<int, int> *kmerTable;
	std::vector<int> *forwardList;
	std::vector<int> *backwardList;
	std::string csvFileName;

	/**
	 * Methods
	 */
	void medianSmooth();
	int findMedian(int, int);
	void score();

public:
	ScorerTr(std::string &, int, int, int);
	virtual ~ScorerTr();
	vector<int> *getForwardScores();
	vector<int> *getBackwardScores();

	void scoresFormat(int, int);

	void printForwardScores(std::string fileName);
	void printBackwardScores(std::string fileName);
};

#endif /* SCORERTR_H_ */
