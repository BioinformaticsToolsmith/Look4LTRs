/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * ModulePipeline
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

#include "ScorerTr.h"
#include "Merger.h"
#include "Detector.h"
#include "Matcher.h"
#include "PostProcess.h"
#include "Filter.h"
#include "LtrParameters.h"
#include "Element.h"
#include "RT.h"
#include "DirectedGraph.h"
#include "DeepNesting.h"

#include "../red/Red.h"


#include <string>
#include <filesystem>
#include <iostream>


class ModulePipeline
{
private:
    // Variables
    Red &red;
    Merger *forwardMerger;
    Merger *backwardMerger;
    std::vector<Element> fElement;
    std::vector<Element> bElement;

    std::vector<Stretch*> heapStretchVec;

    int dbID;

    std::string elePath;

    Matcher *mat;

    DeepNesting *dp;

    // Methods
    void writeToDBHelper(std::string filePath, std::vector<Element> & eleVec, int fastaID, std::string chromName, std::string *chromosome);

    void writeElementsHelper(std::string filePath, std::vector<Element> & eleVec);

    std::vector<Element> readElementsHelper(std::string filePath);

public:
    
    // Constructor

    ModulePipeline(Red &_red);
    ~ModulePipeline();

    // Getter and Setters
    std::vector<RT*>* getRtVec();
    std::vector<RT*>* getComplexVec();
    std::pair<int, int> getFamilyRegion(RT* rt);
    DirectedGraph<Element>* getFamilyGraph(RT *rt);

    // Methods
    void buildElements(std::string *chromosome);
    void writeToDB(std::string filePath,  int fastaID, std::string chromName, std::string *chromosome);
    void writeElements(std::string filePath);
    void clearElementVec();

    void readElements();
    void removeElements();
    void matchElements(IdentityCalculator<int32_t> &icStandard, IdentityCalculator<int32_t> &icRecent, std::string *chromosome);
    void findDeepNests(IdentityCalculator<int32_t> &icStandard, IdentityCalculator<int32_t> &icRecent, std::string *chromosome);
    void process(IdentityCalculator<int32_t> &icStandard, std::string *chromosome, bool nestOnly = false);
    void filter(IdentityCalculator<int32_t> &icStandard, std::string *chromosome);
};