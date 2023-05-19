/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * CaseRecent
 * 
 *  Created on: X X, 20XX
 *      Author: Anthony B. Garza.
 * Reviewer:
 *   Purpose: 
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
#include "RTSolo.h"
#include "CaseSolo.h"
#include "LtrParameters.h"
#include "../IdentityCalculator.h"
#include "../red/Red.h"


#include <iostream>

class CaseRecent : public CaseMatcher {
    public:
        CaseRecent(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq);

        void apply(DirectedGraph<Element>& graph, std::vector<Element*>& forwardVec, std::vector<Element*>& backwardVec, int graphIndex);
    
        bool checkOuterRt(Element *leftLtr, Element *nestedLeftLtr, Element *nestedRightLtr, Element *rightLtr);

        Element * findOuterLtr(Element *e, DirectedGraph<Element> &graph);


};