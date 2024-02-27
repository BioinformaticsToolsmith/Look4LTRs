/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Test.cpp
 *
 *  Created on: Sep 21, 2022
 *      Author: Anthony B. Garza.
 *      Edited by Hani Z. Girgis
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

#pragma once

#include "LtrParameters.h"

#include "Stretch.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <assert.h>

// For testing. See if you need it in the distribution version.
#include <fstream>
#include <tuple>
#include <string>
#include <fstream>

class Merger {
public:

	// Constructor
	Merger(std::vector<int> *scores, int minStretch, int maxGap, int simMargin,
			int interruptMargin, bool isForward);

private:
	/**
	 * Variables
	 */

	// Minimum number of continuous non-zero scores to be marked as Keep
	int minStretch;
	// Maximum gap allowed between merging candidates
	int maxGap;
	// Maximum difference in scores between merging candidates
	int simMargin;
	// Maximum difference in scores between neighboring scores to be considered 'interruptive'
	int interruptMargin;
	// Pointer to scores
	std::vector<int> *scores;
	// Are the scores pointing forward?
	bool isForward;


	/**
	 * Methods
	 */

	/**
	 * Find the median height of a stretch without 0's
	 * start: inclusive
	 * end: exclusive
	 */
	int findMedianScore(int start, int end);



	// Merge a vector of Stretchs depending on criteria
	std::vector<Stretch> merge();

	// Merges Stretchs forward and backwards
	void mergeTwoDirections();

	// remove a certain type of Stretch based on its markType and set the scores of the Stretchs region to 0
	std::vector<Stretch> removeStretches(int mark);

	std::vector<Stretch> stretchVec;  // vector of stretches to later merge

public:

	// Getter
	std::vector<Stretch>* getStretchVec();

	// public attributes

	// for print testing
	void printStretches();

	// prints stretches as scores to file
	void printScores(std::string path);

};
