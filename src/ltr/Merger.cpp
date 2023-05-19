/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Test.cpp
 *
 *  Created on: Sep 21, 2022
 *      Author: Anthony B. Garza.
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
 *
 *
 * To Do:
 * + Clean up
 * + Change at() to []
 */

#include "Merger.h"

Merger::Merger(std::vector<int> *scores, int minStretch, int maxGap,
		int simMargin, int interruptMargin, bool isForward) {

	// Pre-conditions
	assert(scores->size() > 0);
	assert(minStretch > 0);
	assert(maxGap >= 0);
	assert(simMargin >= 0);
	assert(interruptMargin >= 0);

	this->scores = scores;
	this->minStretch = minStretch;
	this->maxGap = maxGap;
	this->simMargin = simMargin;
	this->interruptMargin = interruptMargin;
	this->isForward = isForward;

	// Filling stretches and assigning Keep or Delete
	int prevScore = scores->at(0); // previous score in scores
	int stretchStart = 0; // Start of continuous stretch of scores in scores
	int mark = -1; // Keep or Delete, used from Stretch class; Stretch::K, Stretch::D

	for (int i = 1; i < scores->size(); i++) {
		if (scores->at(i) != prevScore) { // Has the score changed?
			if (prevScore != 0) { // Is this a Stretch, i.e., non-zero stretch of scores?
				mark = i - stretchStart >= minStretch ? Stretch::K : Stretch::D;
				Stretch stretch{stretchStart, i, mark, isForward};
				stretchVec.push_back(stretch);
			}

			prevScore = scores->at(i);
			stretchStart = i;
		}
	}

	// Handling if last stretch reaches to the last score
	if (prevScore != 0) {
		mark = scores->size() - stretchStart >= minStretch ?
				Stretch::K : Stretch::D;
		stretchVec.push_back(
				Stretch { stretchStart, static_cast<int>(scores->size()), mark, isForward });
	}

	if (stretchVec.size() > 0) {

		// First forward and backward merge pass
		mergeTwoDirections();

		// Marking interrupts
		for (size_t i = 1; i < stretchVec.size() - 1; i++) {
			if (stretchVec[i].getMark() == Stretch::D) {

				int medianOfPrev = findMedianScore(stretchVec[i - 1].getStart(),
						stretchVec[i-1].getEnd());
				int medianOfCurr = findMedianScore(stretchVec[i].getStart(),
						stretchVec[i].getEnd());
				int medianOfNext = findMedianScore(stretchVec[i + 1].getStart(),
						stretchVec[i + 1].getEnd());

				if (std::abs(medianOfCurr - medianOfPrev) > interruptMargin
						&& std::abs(medianOfCurr - medianOfNext)
								> interruptMargin) {
					stretchVec[i].setMark(Stretch::I);
				}
			}
		}

		// Removing interrupts
		stretchVec = removeStretches(Stretch::I);

		// Second forward and backward merge pass after removing interrupts
		mergeTwoDirections();

		// Removing all delete Stretches
		stretchVec = removeStretches(Stretch::D);

		// Final forward merge; backward merge is not needed
		stretchVec = merge();

		// Adjusting scores to reflect median heights of stretches and setting height of stretches
		for (int i = 0; i < stretchVec.size(); i++) {
			int start = stretchVec[i].getStart();
			int end = stretchVec[i].getEnd();
			int median = findMedianScore(start, end);
			stretchVec[i].setMedianHeight(median);

			for (int j = start; j < end; j++) {
				scores->at(j) = median;
			}
		}
	}
}


/**
 * Given a vector of scores, a start, and an end location,
 * get the median score not counting scores of 0.
 * Return this median score
 */
int Merger::findMedianScore(int start, int end) {
	// Pre-conditions
	assert(start < end);
	assert(start >= 0);
	assert(end >= 0);
	assert(end <= scores->size());

	// value to return
	int median;

	// creating vector with scores from start to end without 0's
	std::vector<int> stretch;
	for (int i = start; i < end; i++) {
		if (scores->at(i) != 0) {
			stretch.push_back(scores->at(i));
		}
	}
	int stretchSize = stretch.size();

	if (stretchSize == 0) {
		std::cerr << "Region provided in findMedians has no non-zero scores!"
				<< std::endl;
		throw std::exception();
	}

	std::sort(stretch.begin(), stretch.end());
	int middle = stretchSize / 2;
	if (stretchSize % 2 == 1) {
		median = stretch.at(middle);
	}
	else {
		median = int((stretch.at(middle) + stretch.at(middle - 1)) / 2);
	}


	// Post-condition
	assert(median > 0);

	return median;
}

/**
 * Merge stretches depending on the gap between consecutive stretches
 * and the similarity of their scores/heights.
 * Driver: Mr. Anthony B. Garza, Reviewer: Dr. Hani Z. Girgis
 */
std::vector<Stretch> Merger::merge() {
	std::vector<Stretch> mergedVec; // merged Stretches to return
	int i = 0;
	int j = 0;

	int size = static_cast<int>(stretchVec.size());
	while (i < size) {
		Stretch currStretch = stretchVec[i];

		// This is a keep stretch
		if (currStretch.getMark() == Stretch::K) {
			j = i + 1;

			while (j < size) {
				// Checking gap and similar score criteria
				Stretch nextStretch = stretchVec[j];

				if (currStretch.calculateGap(nextStretch) > maxGap
						|| std::abs(
								findMedianScore(currStretch.getStart(),
										currStretch.getEnd())
										- findMedianScore(
												nextStretch.getStart(),
												nextStretch.getEnd()))
								> simMargin) {
					break;
				}

				// If we get here, then we can merge
				currStretch = currStretch.merge(nextStretch, Stretch::K);
				j++;
			}

			mergedVec.push_back(currStretch);
			// Move to the next non-merged Stretch
			i = j;
		}
		// This is a delete stretch
		else {
			mergedVec.push_back(currStretch);

			// Move to the next Stretch if the current one is not a keep Stretch
			i++;
		}
	}

	mergedVec.shrink_to_fit();
	return mergedVec;
}

/**
 * Merge forward, reverse stretches, merge backwards, then reverse stretches back to normal
 */
void Merger::mergeTwoDirections() {
	stretchVec = merge();
	std::reverse(stretchVec.begin(), stretchVec.end());

	stretchVec = merge();
	std::reverse(stretchVec.begin(), stretchVec.end());
}

/**
 * Given the type of stretch to delete (mark), iterate through stretches and if the 
 * type of stretch matches the one to delete then flatten the scores of that stretch's 
 * region to 0 and add the stretch to a vector. Return the vector without the marked stretches  
 */
std::vector<Stretch> Merger::removeStretches(int mark) {
	std::vector<Stretch> removedVector;
	for (auto s : stretchVec) {
		if (s.getMark() == mark) {
			for (int j = s.getStart(); j < s.getEnd(); j++) {
				scores->at(j) = 0;
			}
		} else {
			removedVector.push_back(s);
		}
	}
	removedVector.shrink_to_fit();
	return removedVector;
}

std::vector<Stretch>* Merger::getStretchVec() {
	return &stretchVec; 
}



void Merger::printStretches() {
	std::cout << "Stretches" << std::endl << "------------------------------" << std::endl;
	for (auto s: stretchVec) {
		std::cout << s << " " << findMedianScore(s.getStart(), s.getEnd()) << std::endl;
	}
}

void Merger::printScores(std::string path) {
	std::ofstream file;
	file.open(path);
	int prev_end = 0;
	for (auto s: stretchVec) {

		for (int i = prev_end; i < s.getStart(); i++) {
			file << 0 << std::endl;
		}
		for (int i = 0; i < s.getSize(); i++) {
			file << s.getMedianHeight() << std::endl;
		}
		prev_end = s.getEnd();
	}
	for (int i = prev_end; i < scores->size(); i++) {
		file << 0 << std::endl;
	}

	file.close();
}