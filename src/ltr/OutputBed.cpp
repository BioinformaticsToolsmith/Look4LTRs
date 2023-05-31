/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Bed
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

#include "OutputBed.h"

OutputBed::OutputBed(std::string outputDir, std::string fastaPath) : Output(outputDir, fastaPath)
{
    this->outputPath += ".bed";
    headerVec.insert(headerVec.begin(), {"chrom", "Start", "End"});
}

void OutputBed::write(std::string chrom, std::vector<RT*> &rtVec) {
    if (!isCreated) {
        createFile();
        isCreated = true;
    }

    std::ofstream outFile;
    outFile.open(outputPath, std::ios::app);

    if (outFile.is_open()) {
        if (!rtVec.empty()) {
            outFile << std::endl;

            for (int i = 0; i < rtVec.size() - 1; i++) {
                if (rtVec[i]->hasRightLTR()) {
                    outFile << chrom << "\t" << rtVec[i]->getStart() << "\t" << rtVec[i]->getEnd() << std::endl;
                }
            }
            if (rtVec.back()->hasRightLTR()) {
                outFile << chrom << "\t" << rtVec.back()->getStart() << "\t" << rtVec.back()->getEnd();
            }
        }
        outFile.close();
    }
    else {
        std::cerr << "Could not open " + outputPath + "for writing!" << std::endl;
        throw std::exception();    
    }
}