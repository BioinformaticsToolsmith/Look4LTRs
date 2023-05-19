/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * PostProcess
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

#include "RT.h"
#include "LtrParameters.h"
#include "../IdentityCalculator.h"
#include "LtrUtility.h"

#include <stack>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>

#include <iostream>

#pragma once

class PostProcess
{
private:
    // Variables
    std::vector<RT *> &rtVec;
    std::vector<RT *> &complexVec;
    int max;
    IdentityCalculator<int32_t> &ic;
    std::string *seq;
    // Methods
    void nest();
    void extendIdentity(std::vector<RT *> &vec);
    void extendK(std::vector<RT *> &vec);

    void extendLTR(RT *rtPtr);

public:
    
    // Constructor

    PostProcess(std::vector<RT *> &rtVec, std::vector<RT *> &complexVec, int _max, IdentityCalculator<int32_t> &_ic, std::string *_seq);
    ~PostProcess();

    // Getter and Setters

    // Methods
    void apply(bool nestOnly = false);
};