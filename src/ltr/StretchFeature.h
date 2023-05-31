/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * StretchFeature.h
 * 
 *  Created on: Oct 5, 2022
 *      Author: Anthony B. Garza
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

#include "../Matrix.h"
#include "Stretch.h"
#include "LtrParameters.h"
#include "LtrUtility.h"
#include "../red/Red.h"
#include "../utility/ILocation.h"

#include "assert.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>

class StretchFeature {
public:

	/**
	 * Constructors
	 */

    // This constructor takes the means and the stds of features
	StretchFeature(std::vector<double> &meanVec, std::vector<double> &stdVed);


	/**
	 * Methods
	 */

	// builds a matrix with features from the stretches
	Matrix buildMatrix(const std::vector<Stretch> &stretchVec, Red &red, std::string &seq, std::vector<int> &scoreVec, const std::vector<utility::ILocation *> *locationVec);

	// Scale features using standard scaler and return new Matrix
	Matrix scale(Matrix &featureMatrix);


private:
	/**
	 * Variables
	 */
	std::vector<double> &meanVec;
	std::vector<double> &stdVec;
};
