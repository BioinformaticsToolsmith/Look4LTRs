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

#include "DeepNesting.h"

DeepNesting::DeepNesting(std::vector<RT*>* _rtVecPtr, const std::string *_chrom, ModulePipeline &_mp, Red &_red, IdentityCalculator<int32_t> &_icStandard, IdentityCalculator<int32_t> &_icRecent) : rtVecPtr(_rtVecPtr), chrom(_chrom), mp(_mp), red(_red), icStandard(_icStandard), icRecent(_icRecent)
{

}

DeepNesting::~DeepNesting()
{
    for (auto &mp : mpVec) {
        delete mp;
        mp = nullptr;
    }
}

void DeepNesting::findRegion() {
    for (auto rt : *rtVecPtr) {
        // Searching for the recently nested inside element
        if (rt->getCaseType() == "RecentlyNestedInner") {

            // Get sequence of the graph that these recently nested elements are in
            std::pair<int, int> graphRegion = mp.getFamilyRegion(rt); // start and end of graph
            std::string graphSeq = chrom->substr(graphRegion.first, graphRegion.second - graphRegion.first);

            // Where the nested element starts in the graph
            int nestStart = rt->getStart() - graphRegion.first;
            // How long is the nested element?
            int length = rt->getSize();

            auto nestModuleVec = findDeep(graphSeq, nestStart, length, 2);
            mpVec.insert(mpVec.end(), nestModuleVec.begin(), nestModuleVec.end());

            if (nestModuleVec.size() > 0) {
                // We found a recently nested element!
                auto recentVec = nestModuleVec.at(0)->getRtVec();
                int graphGroup = rt->getGraphGroup();
                auto graph = mp.getFamilyGraph(rt);
            

                // Removing duplicates
                for (auto& r : *recentVec) {
                    r->push(graphRegion.first);
                    bool foundR = false;
                    for (auto &elePtr : graph->getValueVec()) {
                        if (r->hasRightLTR() && (r->getLeftLTR()->calcOverlap(*elePtr) > 0 || r->getRightLTR()->calcOverlap(*elePtr) > 0)) {
                            foundR = true;
                            break;
                        }
                    }

                    

                    if (!foundR) {
                        delete r;
                        r = nullptr;
                    }
                    else {
                        // delete r;
                        // r = nullptr;
                        r->setGraphGroup(graphGroup);
                        recentNestVec.push_back(r);
                    }
                }

            }
        }
    }

    rtVecPtr->insert(rtVecPtr->end(), recentNestVec.begin(), recentNestVec.end());

}

std::vector<ModulePipeline*> DeepNesting::findDeep(std::string &graphSeq, int removeStart, int removeLength, int level) {
    std::string cutSeq = graphSeq.substr(0, removeStart) + graphSeq.substr(removeStart + removeLength);
    ModulePipeline *mp = new ModulePipeline{red};
    std::vector<ModulePipeline*> mpVec{mp};

    mp->buildElements(&cutSeq);
    mp->matchElements(icStandard, icRecent, &cutSeq);
    mp->findRTs();
    mp->process(icStandard, &cutSeq, true);

    // Removing complex
    for (auto cpx : *mp->getComplexVec()) {
        delete cpx;
        cpx = nullptr;
    }

    auto rtVec = mp->getRtVec();
    // for (auto &rt : *rtVec) {
    //     delete rt;
    // }
    std::vector<RT*> r;

    for (auto &rt : *rtVec) {
        // Searching for the recently nested found inside
        if (rt->getCaseType() == "RecentlyNestedInner") {
            #pragma omp critical
            {
                std::cout << "Found a nest at level " << level << std::endl;
            }
            auto nestModulePipeline = findDeep(cutSeq, rt->getStart(), rt->getSize(), level + 1);
            auto nestVec = nestModulePipeline.front()->getRtVec();
            r.insert(r.end(), nestVec->begin(), nestVec->end());
            mpVec.insert(mpVec.end(), nestModulePipeline.begin(), nestModulePipeline.end());

        }

    }
    rtVec->insert(rtVec->end(), r.begin(), r.end());
    for (auto &rt : *rtVec) {
        rt->expand(removeStart, removeLength);
    }

    rtVec->shrink_to_fit();

    LtrUtility::removeNests(*rtVec);
    return mpVec;
}