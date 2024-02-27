/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * DeepNesting
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

#include "ModulePipeline.fwd.h"
#include "DeepNesting.fwd.h"

#include "ModulePipeline.h"
#include "RT.h"
#include "../red/Red.h"
#include "../IdentityCalculator.h"

#include <vector>
#include <string>


class DeepNesting
{
private:
    // Variables
    std::vector<RT *> *rtVecPtr;
    const std::string *chrom;
    ModulePipeline &mp;
    std::vector<ModulePipeline*> mpVec;
    std::vector<RT*> recentNestVec;
    Red &red;
    IdentityCalculator<int32_t> &icStandard;
    IdentityCalculator<int32_t> &icRecent;

    // Methods

public:
    
    // Constructor

    DeepNesting(std::vector<RT*>* _rtVecPtr, const std::string *_chrom, ModulePipeline &_mp, Red &_red, IdentityCalculator<int32_t> &_icStandard, IdentityCalculator<int32_t> &_icRecent);
    ~DeepNesting();

    // Methods
    void findRegion();
    std::vector<ModulePipeline*> findDeep(std::string &graphSeq, int removeStart, int removeLength, int level);
};