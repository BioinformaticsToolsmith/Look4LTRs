/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Output
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

Output::Output(std::string outputDir, std::string fastaPath)
{
    if (!std::filesystem::exists(outputDir)) {
        std::filesystem::create_directory(outputDir);
    }
    this->outputPath = outputDir + "/" + std::filesystem::path(fastaPath).filename().stem().string();

}

void Output::createFile() {
    std::ofstream outFile;
    outFile.open(outputPath, std::ios::trunc);

    if (outFile.is_open()) {
        std::copy(headerVec.begin(), headerVec.end(), std::ostream_iterator<std::string>(outFile, "\t"));
        outFile.close();
    }
    else {
        std::cerr << "Could not open " + outputPath + "for writing!" << std::endl;
        throw std::exception();
    }
}

