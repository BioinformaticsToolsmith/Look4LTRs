/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * SGD.cpp
 * 
 *  Created on: Oct 5, 2022
 *      Author: Anthony B. Garza
 *      Reviewer: Hani Z. Girgis, PhD
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

#include "SGD.h"

SGD::SGD(std::vector<double> weightVec) {
	// Pre-conditions
	assert(weightVec.size() > 0);  // Number of features must be greater than 0!

	this->weightVec = weightVec;
	this->featureCount = weightVec.size() - 1;

}


std::vector<int> SGD::predict(Matrix &featureMatrix) {

	int rowCount = featureMatrix.getNumRow();
	std::vector<int> predictionVec(rowCount, 0);

	// Looping through each instance
	for (int row = 0; row < rowCount; row++) {
		// Will contain final prediction of current instance
		double pred = weightVec[0];

		// For each feature, multiply its weight and add to pred
		// Note that weightVec starts with the bias; add 1 to get the feature weight
		for (int col = 0; col < featureCount; col++) {
			pred += featureMatrix(row, col) * weightVec[col + 1];
		}
		predictionVec[row] = pred <= 0 ? 0 : 1;
	}

	return predictionVec;
}

std::vector<double> SGD::getWeightVec() {
	return weightVec;
}

void SGD::setWeightVec(std::vector<double>& weightVec) {
	this->weightVec = weightVec;
	this->featureCount = weightVec.size() - 1;
}
