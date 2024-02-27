/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * LtrParameters.cpp
 * 
 *  Created on: Oct 12, 2022
 *      Author: Hani Z. Girgis, PhD and Anthony B. Garza.
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


#include "LtrParameters.h"

std::vector<double> LtrParameters::WEIGHT_VECTOR( { 
			-1.28120866,
			-0.0010650548016567808,
			-0.0012076114812354134,
			-2.663910196609594,
			-0.0010688756484007384,
			-0.0005200200424128672,
			0.0069863391290507434,
			-0.001116804159415559,
			-0.001753367206601686,
			0.0007393122704765399,
			-0.0031334346794077757} );

std::vector<double> LtrParameters::MEAN_VECTOR( {
			991.7131036238649,
			990.9998302639395,
			167815.26483068828,
			3764.6503012815074,
			457.11743613680727,
			660.4753328022987,
			681.85360944034,
			681.8938962513797,
			168.40618045777393} );

std::vector<double> LtrParameters::STD_VECTOR( { 
			1059.009702424861,
			1058.4686450779059,
			346956.290142164,
			4535.429654736708,
			613.4340387344209,
			4100.770401833956,
			2960.003971700431,
			2960.0085150338346,
			4444.491631354879} );

double LtrParameters::MIN_PERC = 0.3266381448873077;

double LtrParameters::MIN_WEIGHT = 0.26802415;



void LtrParameters::updateParameters(const std::string& filePath, bool hasMatcherValues) {
	std::ifstream file(filePath);
	std::string line;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string token;
		iss >> token; // Read the prefix

		std::vector<double> *targetVector = nullptr;

		if (token == "mean:") {
			targetVector = &MEAN_VECTOR;
		} else if (token == "std:") {
			targetVector = &STD_VECTOR;
		} else if (token == "weights:") {
			targetVector = &WEIGHT_VECTOR;
		} else if (token == "intercept:") {
			// Since intercept is part of WEIGHT_VECTOR, we handle it separately
			WEIGHT_VECTOR.clear();
			double intercept;
			iss >> intercept; // Read the intercept value
			WEIGHT_VECTOR.push_back(intercept);
			continue;
		} else if (hasMatcherValues) {
			if (token == "red") {
				iss >> MIN_PERC;
			} else if (token == "weight") {
				iss >> MIN_WEIGHT;
			}
		}

		if (targetVector != nullptr) {
			// Clear existing values
			targetVector->clear();
			double value;
			while (iss >> value) {
				targetVector->push_back(value);
			}
		}
	}

	file.close();
	
}
