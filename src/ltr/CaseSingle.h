/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * SingleCase
 * 
 *  Created on: Oct 27, 2022
 *      Author: Anthony B. Garza.
 *    Reviewer: Hani Z. Girgis
 *     Purpose: Find (or not) a single complete LTR RT 
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

#include "CaseMatcher.h"
#include "RTComplete.h"
#include "../red/Red.h"
#include "../IdentityCalculator.h"
#include "LtrUtility.h"
#include "unordered_map"

// FOR TESTING @@@@@@@@@@@@@@@@@@@@@@ DELETE @@@@@@@@@@@@@@@@@@@@@@@
#include <iostream>

class CaseSingle: public CaseMatcher
{
    public:
        CaseSingle(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq);

        void apply(DirectedGraph<Element>& graph, std::vector<Element*>& forwardVec, std::vector<Element*>& backwardVec, int graphIndex);

        bool checkSingle(Element* left, Element* right);
        bool checkSingle(std::vector<int> &leftLtrScoreVec, std::vector<int> &interiorScoreVec, std::vector<int> &rightLtrScoreVec);
        void filterSequential();



};

