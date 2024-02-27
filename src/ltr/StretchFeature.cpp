/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * StretchFeature.cpp
 * 
 *  Created on: Oct 5, 2022
 *      Author: Anthony B. Garza
 *    Reviewer: Hani Z. Girgis
 * 
 *  Purpose: Contains the means and standard deviations of each feature used for merging consecutive stretches
 *           Is capable of building a matrix of these features as well as scaling a matrix of features
 *           according to the normalization of (x - mean) / standard_deviation
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

#include "StretchFeature.h"

//StretchFeatures::StretchFeatures() {
//	isMeans = false;
//	isStds = false;
//
//	// Giving default value of 0 to means and stds
//	means.fill(0.);
//	stds.fill(0.);
//}
//

StretchFeature::StretchFeature(std::vector<double> &meanVec,
		std::vector<double> &stdVec) :
		meanVec(meanVec), stdVec(stdVec) {
	assert(meanVec.size() > 0);
	assert(stdVec.size() > 0);
	assert(meanVec.size() == stdVec.size());

}

Matrix StretchFeature::buildMatrix(const std::vector<Stretch> &stretchVec, Red &red, std::string &seq, std::vector<int> &scoreVec, const std::vector<utility::ILocation *> *locationVec) {
	// The matrix size is the number of stretches minus 1 (because features are
	// extracted from consecutive pairs) times the features.
	
	/**
	 * There are 10 features:
	 * First Stretch Size     : size of the first stretch
	 * Second Stretch Size    : size of the second stretch
	 * Gap Size               : size of the gap between the stretches
	 * Height Difference      : difference of stretches' heights
	 * Median Score Difference: difference of median (does not count zeros) Red scores of both stretches
	 * Mean Score Difference  : difference of mean (involves zero in calculations) Red scores of both stretches
	 * First Mean Score       : mean Red Score of the first stretch
	 * Second Mean Score      : mean Red Score of the second stretch
	 * Gap Mean Score         : mean Red Score of the gap between the stretches
	 * Repeat overlap         : 1 if both stretches overlap the same repeat detected by Red
	*/

	Matrix r;

	if (!stretchVec.empty()) {

		r = Matrix(stretchVec.size() - 1, 10);

		std::vector<int> overlapVec(stretchVec.size() - 1, -1);
		int i = 0; 
		int j = 0;
		int label = -1;
		while (i < stretchVec.size() - 1) {
			if (j < locationVec->size()) {
				if (stretchVec[i].getEnd() < locationVec->at(j)->getStart()) {
					label = 0;
				}
				else if (stretchVec[i + 1].getStart() > locationVec->at(j)->getEnd() + 1){
					j++;
					continue;
				}
				else if (stretchVec[i].isPartOf(locationVec->at(j)->getStart(), locationVec->at(j)->getEnd() + 1, LtrParameters::MIN_STRETCH_OVERLAP) && 
						stretchVec[i + 1].isPartOf(locationVec->at(j)->getStart(), locationVec->at(j)->getEnd() + 1, LtrParameters::MIN_STRETCH_OVERLAP)){
					label = 1;
				}
				else {
					label = 0;
				}
			}
			else {
				label = 0;
			}

			overlapVec[i] = label;
			i++;

		}

		

		for (int i = 0; i < stretchVec.size() - 1; i++) {
			r(i, 0) = stretchVec[i].getSize();
			r(i, 1) = stretchVec[i + 1].getSize();
			r(i, 2) = stretchVec[i].calculateGap(stretchVec[i + 1]);
			r(i, 3) = std::abs(
					stretchVec[i].getMedianHeight()
							- stretchVec[i + 1].getMedianHeight());

			std::vector<int> firstStretch(scoreVec.begin() + stretchVec[i].getStart(), scoreVec.begin() + stretchVec[i].getEnd());
			std::vector<int> secondStretch(scoreVec.begin() + stretchVec[i + 1].getStart(), scoreVec.begin() + stretchVec[i + 1].getEnd());
			std::vector<int> gapStretch(scoreVec.begin() + stretchVec[i].getEnd(), scoreVec.begin() + stretchVec[i + 1].getStart());

			double firstMean =  LtrUtility::calcMean(firstStretch);
			double secondMean =  LtrUtility::calcMean(secondStretch);
			double gapMean =  gapStretch.size() > 0? LtrUtility::calcMean(gapStretch) : 0.0;

			r(i, 4) =  std::abs(LtrUtility::calcMedian(firstStretch) - LtrUtility::calcMedian(secondStretch));
			r(i, 5) =  std::abs(firstMean - secondMean);
			r(i, 6) = firstMean;
			r(i, 7) = secondMean;
			r(i, 8) = gapMean;
			r(i, 9) = overlapVec[i];
		}

	}

	return r;
}


Matrix StretchFeature::scale(Matrix &m) {
	int rowCount = m.getNumRow();
	int colCount = m.getNumCol();

	Matrix r(rowCount, colCount);  // result

	// All features but the Repeat Overlap
	for (int col = 0; col < colCount - 1; col++) {
		for (int row = 0; row < rowCount; row++) {
			r(row,col) = (m(row,col) - meanVec[col]) / stdVec[col];
		}
	}

	return r;
}

