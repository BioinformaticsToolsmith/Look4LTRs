/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * CaseRecentComplex
 * 
 *  Created on: Dec 15, 2022
 *      Author: Anthony B. Garza.
 * Reviewer:
 *   Purpose: If all other cases fail, assume that all elements are solo LTRs
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

#include "Element.h"
#include "CaseMatcher.h"
#include "RTComplex.h"
#include "DirectedGraph.h"
#include "../red/Red.h"
#include "../IdentityCalculator.h"
#include "CaseSingle.h"

#include <unordered_set>

// @@@@@@@@@@@@@@@@@@@ FOR TESTING: DELETE AFTER @@@@@@@@@@@@@@@@@@@@@@@
#include <iostream>


class CaseRecentComplex : public CaseMatcher {
    public:
        CaseRecentComplex(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq);

        void apply(DirectedGraph<Element> &graph, std::vector<Element *> &forwardVec, std::vector<Element *> &backwardVec, int graphIndex);
};

