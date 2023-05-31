/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * CaseRecentComplex
 * 
 *  Created on: Dec 15, 2022
 *      Author: Anthony B. Garza.
 *    Reviewer:
 *     Purpose: If all other cases fail, assume that all elements are solo LTRs
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

#include "CaseRecentComplex.h"

CaseRecentComplex::CaseRecentComplex(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq) : CaseMatcher(_ic, _icRecent, _red, _seq){
    name = "RecentComplex";
    rank = -1;
}

void CaseRecentComplex::apply(DirectedGraph<Element> &graph, std::vector<Element *> &forwardVec, std::vector<Element *> &backwardVec, int graphIndex)
{

    /**
     * 
     *
     * > > > ....
     *
     * O  O oO O   O O 
     *  \ |\ |\
     *   \| \| \
     *    O  O ...
     *
     * If the areas between the LTR candidates are mostly repetitive as well as similar to each other, this area is extremely complex
     * and should be noted as such. Another case may fire off in this region, but for all intents and purposes, it is just the best guess.
     */

    CaseSingle cs{ic, icRecent, red, seq};


    // Ensures this is more complicated than a single case, at the very least.
    if (forwardVec.size() > 1 && backwardVec.size() > 1) {
        int identitySum = 0;
        std::vector<std::string> interiorVec;

        for (auto fElePtr : forwardVec) {
            auto conVec = getDiagonals(fElePtr, graph);
            Element* bElePtr = nullptr;
            for (auto ele : conVec) {
                if (ele->getStart() > fElePtr->getEnd() && graph.isBidirectional(*fElePtr, *ele)) {
                    bElePtr = ele;
                    break;
                }
            }

            if (bElePtr != nullptr && cs.checkSingle(fElePtr, bElePtr)) {
                interiorVec.push_back(substrInterior(fElePtr, bElePtr, seq));
            }
        }

        // All vs All
        if (interiorVec.size() > 1){
            for (int i = 0; i < interiorVec.size() - 1; i++) {
                for (int j = i + 1; j < interiorVec.size(); j++) {
                    if (areSeqSame(interiorVec.at(i), interiorVec.at(j))) {
                        identitySum++;
                    }
                }
            }
        }

        interiorVec.clear();
        if (identitySum > 1) {

            std::unordered_set<Element *> elementSet;
            std::vector<Element*> eleVec;
            for (auto elePtr : graph.getValueVec()) {
                // If looking at a forward element with a vertical connection to a backwards, we don't want to revisit the backwards.
                if (elementSet.count(elePtr) == 0) {
                    auto vert = getVertical(elePtr, graph);
                    Element *ele = nullptr;

                    // Merge the element and its vertical if it has one.
                    if (vert != nullptr) {
                        ele = new Element{*elePtr, *vert};
                        elementSet.insert(vert);
                    }
                    else {
                        ele = elePtr;
                    }

                    eleVec.push_back(ele);
                }
            }
            std::sort(eleVec.begin(), eleVec.end(), [](Element *e1, Element *e2){return e1->getStart() < e2->getStart();});
            rtVec.push_back(new RTComplex(eleVec, name, rank, graphIndex));
        }
    
    }


}