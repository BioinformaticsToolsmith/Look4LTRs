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

#include "OutputCpx.h"

OutputCpx::OutputCpx(std::string outputDir, std::string fastaPath) : Output(outputDir, fastaPath)
{
    this->outputPath += ".cpx";
    headerVec.insert(headerVec.begin(), {"chrom", "GraphID", "Start", "End", "..."});
}

void OutputCpx::write(std::string chrom, std::vector<RT*> &cpxVec) {
    if (!isCreated) {
        createFile();
        isCreated = true;
    }

    std::ofstream outFile;
    outFile.open(outputPath, std::ios::app);

    if (outFile.is_open()) {
        if (!cpxVec.empty()) {
            outFile << std::endl;

            for (int i = 0; i < cpxVec.size() - 1; i++) {
                outFile << chrom << "\t" << cpxVec.at(i)->getGraphGroup() << "\t" << cpxVec.at(i)->getStart() << "\t" << cpxVec.at(i)->getEnd();
                for (auto ele : cpxVec.at(i)->getLTRVec()) {
                    outFile << "\t" << ele->getStart() << "\t" << ele->getEnd();
                }
                outFile << std::endl;
            }
            outFile << chrom << "\t" << cpxVec.back()->getStart() << "\t" << cpxVec.back()->getEnd();
            for (auto ele : cpxVec.back()->getLTRVec()) {
                outFile << "\t" << ele->getStart() << "\t" << ele->getEnd();
            }
        }
        outFile.close();
    }
    else {
        std::cerr << "Could not open " + outputPath + "for writing!" << std::endl;
        throw std::exception();    
    }
}