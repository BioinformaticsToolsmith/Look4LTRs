/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * OutputCpx
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

#include "Output.h"
#include "RT.h"

class OutputCpx : Output
{
private:
    // Variables

    // Methods

public:
    
    // Constructor

    OutputCpx(std::string outputDir, std::string fastaPath);
    void write(std::string chrom, std::vector<RT*> &cpxVec);

    // Getter and Setters

    // Methods
};