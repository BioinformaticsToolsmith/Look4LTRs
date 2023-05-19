/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Filter
 * 
 *  Created on: 12 6, 2022
 *      Author: Anthony B. Garza.
 * Reviewer:
 *   Purpose: For filtering out LTR RT false positives;
 *            Filters by structural characteristics, 
 *            repetivity in the genome, and clustering
 *        
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

#include "RT.h"
#include "LtrParameters.h"
#include "LtrUtility.h"

#include "../red/Red.h"
#include "LocalAlignment.h"

#include <vector>
#include <stack>
#include <assert.h>
#include <algorithm>
#include <unordered_set>

#include <iostream>

typedef std::vector<std::pair<int, int>> Range;


class Filter
{
private:
    // Variables
    std::vector<RT*>& rtVec;
    Red& red;
    std::string *seq;
    // Methods


    // Filters RTs by the minimum and maximum lengths for the LTRs, interiors, and RT LTRs.
    void filterLength();

    // void filterLengthRatio();


    void filterIdentity();

    // Filters RTs by their repetivity. If not repetitive throughout the genome, they are not TEs.
    void filterRepetivity();

    bool findMITE(std::string &seq);

    // Filters LTRs by searching for MITEs
    void filterMITE();

    // Helper method to find a PPT in a sequence for filterPPT
    std::pair<int, int> findPPT(int border, int size, char replaceChar, char withChar);

    // Filters RTs by their poly purine trail.  It is expected for interiors to be dropped here.
    void filterPPT();

    // marks if a RT has a TSD
    void markTSD();

    // Filter solo LTRs without a graph mate
    void filterLonesomes();

    // Filters RTs by clustering
    void filterCluster();

    // If a RT is to be filtered out, remove all connections to and fro other RTs
    void removeNests(RT *rt);

    std::vector<int> scoreSeq(int s, int e);



public:
    
    // Constructor

    // rtVec should be a vector of RT pointers, i.e., the candidate LTR RTs.
    Filter(std::vector<RT*> &_rtVec, Red &_red, std::string *_seq);

    // Getter and Setters
    std::vector<RT*> getRtVec();

    // Methods
    void apply();
};