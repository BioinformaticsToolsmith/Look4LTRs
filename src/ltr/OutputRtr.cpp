/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * OutputRTR
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

#include "OutputRtr.h"

OutputRtr::OutputRtr(std::string outputDir, std::string fastaPath) : Output(outputDir, fastaPath)
{
    this->outputPath += ".rtr";
    headerVec.insert(headerVec.begin(), {"chrom", "ID", "LeftStart", "LeftEnd", "RightStart", "RightEnd", "NestIn", "NestOut", "RC", "PPTStart", "PPTEnd", "TSDStart", "TSDEnd", "CaseType", "GraphGroup", "LTRIdentity"});
}

std::string OutputRtr::formatRt(std::string& chrom, RT* rt, std::unordered_map<RT*, int> &idTable) {
    std::stringstream s;
    s << chrom << "\t" << idTable[rt] << "\t" << rt->getStart() << "\t";
    // Is a CompleteRT
    if (rt->hasRightLTR()) {
        s << rt->getLeftLTR()->getEnd() << "\t" << rt->getRightLTR()->getStart() << "\t" << rt->getEnd() << "\t";
    }
    else {
        s << rt->getEnd() << "\t" << "NA" << "\t" << "NA" << "\t";
    }

    // Has nested?
    if (rt->hasRightLTR() && rt->hasNest()) {
        s << "{";
        for (auto nestPtr : rt->getNestSet()) {
            s << idTable[nestPtr] << ",";
        }
        s << "}" << "\t";
    }
    else {
        s << "NA" << "\t";
    }

    // Is nested?
    if (rt->hasOuter()) {
        s << idTable[rt->getOuterNest()] << "\t";

    }
    else {
        s << "NA" << "\t";
    }

    if (rt->hasRightLTR()) {
        std::string rc = rt->getIsRC() ? "1" : "0";
        s << rc << "\t";

        s << rt->getPPTStart() << "\t" << rt->getPPTEnd() << "\t";

        if (rt->getIsTSDExist()) {
            s << rt->getLeftTSD().first << "\t" << rt->getRightTSD().second << "\t";
        }
        else {
            s << "NA" << "\t" << "NA" << "\t";
        }
    }
    else {
        s << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t";
    }
    s << rt->getCaseType() << "\t" << rt->getGraphGroup();
    if (rt->hasRightLTR()) {
        s << "\t" << rt->getIdentityScore();
    }
    else {
        s << "\t" << "NA";
    }

    return s.str();
}

void OutputRtr::write(std::string chrom, std::vector<RT*> &rtVec) {
    if (!isCreated) {
        createFile();
        isCreated = true;
    }

    std::ofstream outFile;
    outFile.open(outputPath, std::ios::app);

    if (outFile.is_open()) {

        if (!rtVec.empty()) {
            // Building ID Table; For retaining nested element relationships
            std::unordered_map<RT*, int> idTable;

            for (auto rtPtr : rtVec) {
                idTable[rtPtr] = id;
                id++;
            }

            outFile << std::endl;

            for (int i = 0; i < rtVec.size() - 1; i++) {
                outFile << formatRt(chrom, rtVec[i], idTable);
                outFile << std::endl;
            }
            outFile << formatRt(chrom, rtVec.back(), idTable);
        }
        outFile.close();
    }
    else {
        std::cerr << "Could not open " + outputPath + "for writing!" << std::endl;
        throw std::exception();    
    }
}