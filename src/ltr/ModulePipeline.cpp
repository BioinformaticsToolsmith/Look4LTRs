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

#include "ModulePipeline.h"

ModulePipeline::ModulePipeline(Red &_red) : red(_red)
{
    mat = nullptr;
    forwardMerger = nullptr;
    backwardMerger = nullptr;
    dp = nullptr;
}

ModulePipeline::~ModulePipeline()
{
    if (mat != nullptr)
    {
        delete mat;
        mat = nullptr;
    }    

    if (forwardMerger != nullptr)
    {
        delete forwardMerger;
        forwardMerger = nullptr;
    }
    if (backwardMerger != nullptr)
    {
        delete backwardMerger;
        backwardMerger = nullptr;
    }
    for (auto &stretch : heapStretchVec)
    {
        delete stretch;
        stretch = nullptr;
    }
    if (dp != nullptr)
    {
        delete dp;
        dp = nullptr;
    }
}

void ModulePipeline::buildElements(std::string *chromosome)
{
    // Score the sequence
    // std::cout << "Scoring the sequence..." << std::endl;
    ScorerTr *st = new ScorerTr(*chromosome, LtrParameters::K, LtrParameters::MIN_DISTANCE, LtrParameters::MAX_DISTANCE);


    // Merge the scores
    // std::cout << "Merging the scores forward..." << std::endl;
    forwardMerger = new Merger(st->getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP, LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);
    // std::cout << "Size of stretches: " << forwardMerger->getStretchVec()->size() << std::endl;

    // std::cout << "Merging the scores backward..." << std::endl;
    backwardMerger = new Merger(st->getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP, LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);



    // Detect the elements
    // std::cout << "Detecting the elements..." << std::endl;
    Detector dt{red, *chromosome};
    fElement = dt.apply(*forwardMerger->getStretchVec());
    bElement = dt.apply(*backwardMerger->getStretchVec());
    // std::cout << "Done detecting the elements." << std::endl;

}

void ModulePipeline::writeToDBHelper(std::string filePath, std::vector<Element> & eleVec, int fastaID, std::string chromName, std::string *chromosome) {
    std::fstream dbFile{filePath, ios_base::app};
    for (auto &ele : eleVec) {
        if (ele.getSize() >= LtrParameters::MIN_LTR) {
            dbFile << ">" << fastaID << "_" << chromName << "_" << dbID << std::endl;
            dbFile << chromosome->substr(ele.getStart(), ele.getEnd() - ele.getStart()) << std::endl;

            dbID += 1;
        }
    }
    dbFile.close();
}

void ModulePipeline::writeToDB(std::string filePath, int fastaID, std::string chromName, std::string *chromosome) {
    writeToDBHelper(filePath, fElement, fastaID, chromName, chromosome);
    writeToDBHelper(filePath, bElement, fastaID, chromName, chromosome);
}

void ModulePipeline::writeElementsHelper(std::string filePath, std::vector<Element> & eleVec) {
    std::ofstream eleFile{filePath};
    for (auto &ele : eleVec) {
        if (ele.getSize() >= LtrParameters::MIN_LTR) {
            eleFile << ele << std::endl;
        }
    }
    eleFile.close();
}

void ModulePipeline::writeElements(std::string filePath) {
    elePath = filePath;
    writeElementsHelper(elePath + ".fele", fElement);
    writeElementsHelper(elePath + ".bele", bElement);
}


void ModulePipeline::clearElementVec() {
    fElement.clear();
    bElement.clear();
    delete forwardMerger;
    delete backwardMerger;
    forwardMerger = nullptr;
    backwardMerger = nullptr;
}

std::vector<Element> ModulePipeline::readElementsHelper(std::string filePath) {
    std::vector<Element> r;
    std::fstream eleFile{filePath};
    std::string line;
    std::string extra;
    std::string stretchData;

    while (std::getline(eleFile, line)) {
        std::vector<std::tuple<int, int, int>> stretchDataVec;
        bool isForward = true;
        std::stringstream ss{line};
        ss >> extra;
        while (ss >> stretchData) {
            if (stretchData != "-" && stretchData != "+") {
                int medianHeight;
                int start;
                int end;

                std::stringstream ssStretch{stretchData};
                std::string token;

                std::getline(ssStretch, token, ',');
                medianHeight = std::stoi(token.substr(1));

                std::getline(ssStretch, token, ',');
                start = std::stoi(token);

                std::getline(ssStretch, token);
                end = std::stoi(token.substr(0, token.size() - 1));

                stretchDataVec.push_back(std::tuple<int, int, int>{medianHeight, start, end});
                
            }
            else {
                isForward = stretchData == "+"? true:false;
            }
        }
        std::vector<Stretch*> stretchVec;

        for (auto &x : stretchDataVec) {
            Stretch *newStretch = new Stretch{std::get<1>(x), std::get<2>(x), Stretch::K, isForward};
            stretchVec.push_back(newStretch);
            heapStretchVec.push_back(newStretch);
            stretchVec.back()->setMedianHeight(std::get<0>(x));
        }

        Element e{*stretchVec.at(0)};
        int x = 1;
        while (x < stretchVec.size()) {
            e.merge(*stretchVec.at(x));
            x++;
        }

        r.push_back(e);
    }
    return r;
}

void ModulePipeline::readElements() {

    if (!std::filesystem::exists(elePath + ".fele") || !std::filesystem::exists(elePath + ".bele")) {
        std::cout << "Error: Element files not found:" << elePath << std::endl;
        exit(1);
    }

    fElement = readElementsHelper(elePath + ".fele");
    bElement = readElementsHelper(elePath + ".bele");
}

void ModulePipeline::removeElements() {
    std::filesystem::remove(elePath + ".fele");
    std::filesystem::remove(elePath + ".bele");
}

void ModulePipeline::matchElements(IdentityCalculator<int32_t> &icStandard, IdentityCalculator<int32_t> &icRecent, std::string *chromosome ) {

    mat = new Matcher{fElement, bElement, red, icStandard, icRecent, chromosome};
}

void ModulePipeline::findDeepNests(IdentityCalculator<int32_t> &icStandard, IdentityCalculator<int32_t> &icRecent, std::string *chromosome) {
    assert(mat != nullptr);
    process(icStandard, chromosome, true);
    dp = new DeepNesting(mat->getRtVec(), chromosome, *this, red, icStandard, icRecent);
    // // Finding the deep nests
    dp->findRegion();
    auto rtVecPtr = mat->getRtVec();
    // We need to remove the nests so we can renest with the new recently nested elements we found
    LtrUtility::removeNests(*rtVecPtr);
    // Remove any LTR RTs that were ill-formated by the mapping back procedure from recently nested elements
    LtrUtility::removeIllFormat(*rtVecPtr);
    // Ranking the LTR RTs, e.g., if a single and a Recently Nested case were found in the same spot, only report the Recently Nested
    LtrUtility::rankRTs(*rtVecPtr, true);
    // Removing LTR RTs where both LTRs of an RT have large overlap with each other
    LtrUtility::removeDuplicateRTs(*rtVecPtr);
    // Renaming the RecentlyNestedOuter and RecentlyNestedInner to just RecentlyNested. Possible multi-level nesting, after all.
    LtrUtility::renameCaseType(*rtVecPtr, "RecentlyNestedOuter", "RecentlyNested");
    LtrUtility::renameCaseType(*rtVecPtr, "RecentlyNestedInner", "RecentlyNested");

}


void ModulePipeline::process(IdentityCalculator<int32_t> &icStandard, std::string *chromosome, bool nestOnly) {
    PostProcess post{*mat->getRtVec(), *mat->getComplexVec(), chromosome->size(), icStandard, chromosome};
    post.apply(nestOnly);
}

void ModulePipeline::filter(IdentityCalculator<int32_t> &icStandard, std::string *chromosome) {
    Filter filter(*mat->getRtVec(), red, chromosome);
    filter.apply();
    for (auto rt : *mat->getRtVec()) {
        if (rt->hasRightLTR()) {
            std::string leftLTR = chromosome->substr(rt->getLeftLTR()->getStart(), rt->getLeftLTR()->getSize());
            std::string rightLTR = chromosome->substr(rt->getRightLTR()->getStart(), rt->getRightLTR()->getSize());
            rt->setIdentityScore(icStandard.score(&leftLTR, &rightLTR));
        }
    }
}

std::vector<RT*>* ModulePipeline::getRtVec() {
    return mat->getRtVec();
}

std::vector<RT*>* ModulePipeline::getComplexVec() {
    return mat->getComplexVec();
}

std::pair<int, int> ModulePipeline::getFamilyRegion(RT* rt) {
    auto graph = mat->getSubGraph(rt);
    auto vec = graph->getValueVec();
    int start = vec.front()->getStart();
    int end = vec.front()->getEnd();
    for (auto &x : vec) {
        start = std::min(start, x->getStart());
        end = std::max(end, x->getEnd());
    }
    return std::pair<int, int>{start, end};
}

DirectedGraph<Element>* ModulePipeline::getFamilyGraph(RT *rt) {
    return mat->getSubGraph(rt);
}

