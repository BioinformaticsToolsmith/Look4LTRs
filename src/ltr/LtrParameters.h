/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Parameters.h
 * 
 *  Created on: Sep 27, 2022
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
#ifndef SRC_LTR_LTRPARAMETERS_H_
#define SRC_LTR_LTRPARAMETERS_H_

#include <array>
#include <vector>

class LtrParameters {
public:

	static const int K = 13;

	static const int MIN_DISTANCE = 250;

	static const int MAX_DISTANCE = 40000;

	// Minimum size of a stretch to be designated as a 'Keep' stretch
	static const int MIN_STRETCH = 16;

	// Maximum size of a gap allowed between two stretches to be considered for merging
	static const int MAX_GAP = 75;

	// Max difference in height allowed between two stretches to allow for merging
	static const int SIM_MARGIN = 75;

	// Max difference in height between surrounding stretches for the inner stretch to be considered an 'Interrupt' stretch
	static const int INT_MARGIN = 95;

	// The initial value used by the kmer table in the ScorerTr class
	static const int INITIAL_VALUE = -10;

	// The initial value of a score used in the ScorerTr class
	static const int INITIAL_SCORE = 0;

	// Threshold of strong match
	// static constexpr double MATCH_THRESHOLD = 0.5;

	// Threshold for two elements to be overlapping
	static constexpr double OVERLAP_THRESHOLD = 0.1;

	// Minimum percentage of non-zero red scores in an area to be considered a repeat; determined by 2nd percentile
	static constexpr double MIN_PERC = 0.3266381448873077;
	// static constexpr double MIN_PERC = 0.25;


	// Minimum ratio of an LTR's median RED score over its interior's median RED score
	// static constexpr double MIN_LTR_INT_RATIO = 0.5688996168582375;
	// static constexpr double MIN_LTR_INT_RATIO = 0.0;


	static std::vector<double> WEIGHT_VECTOR;

	static std::vector<double> MEAN_VECTOR;

	static std::vector<double> STD_VECTOR;


	// Weight minimum threshold reciprocal for the directed graph
	static constexpr double MIN_WEIGHT = 0.26802415;

	// Minimum overlap needed between stretch and repeat divided by stretch size
	static constexpr double MIN_STRETCH_OVERLAP = 0.50;

	static constexpr double MIN_LENGTH_RATIO = 0.9;

	static constexpr double MIN_IDENTITY_FILTER = 0.5;

	// Min identity score between two sequences
	static constexpr double MIN_IDENTITY = 0.80;

	// Min identity score between two recently nested sequences
	static constexpr double MIN_IDENTITY_RECENT = 0.6;

	// Min size of LTR; inclusive
	static const int MIN_LTR = 80;

	// Max size of LTR; inclusive
	static const int MAX_LTR = 7000;

	// Min size of interior region of RT; inclusive
	static const int MIN_INTERIOR = 70;

	// Max size of interior region of RT; inclusive
	static const int MAX_INTERIOR = 25900;

	// Min size of RT; inclusive
	static const int MIN_RT = 600;

	// Max size of RT; inclusive
	static const int MAX_RT = 34500;

	// MIN MITE Inverted Repeat size
	static const int MIN_MITE_TIR_SIZE = 15;

	// Max MITE Inverted Repeat size
	static const int TIR_SEARCH_RANGE = 30;

	// MAX MITE size
	static const int MAX_MITE_SIZE = 500;

	// How similar do the TIRS of MITES need to be?
	static constexpr double MITE_SIMILARITY = 0.85;

	// Min size of poly purine trail
	static const int MIN_PPT_SIZE = 10;

	// Max size of poly purine trail
	static const int MAX_PPT_SIZE = 100;

	// How far out the PPT can be from the LTR
	static const int MAX_PPT_DISTANCE = 400;

	// MAX TSD Distance
	static const int MAX_TSD_DISTANCE = 20;

	// Min size of TSD
	static const int MIN_TSD_SIZE = 4;

	static constexpr double MAX_N_RATIO = 0.25;
};

#endif /* SRC_LTR_LTRPARAMETERS_H_ */
