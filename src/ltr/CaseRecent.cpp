/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * CaseRecent
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

#include "CaseRecent.h"

CaseRecent::CaseRecent(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq) : CaseMatcher(_ic, _icRecent, _red, _seq){
    name = "RecentlyNested";
    rank = 400;
}

/**
 * [Reviewed]
 */
void CaseRecent::apply(DirectedGraph<Element> &graph, std::vector<Element *> &forwardVec, std::vector<Element *> &backwardVec, int graphIndex) {
    /**
     *
     * >-->-->-->
     * 
     * O  O  O
     *  \ }\ }\
     *   \} \} \
     *    O  O  O
     *  
    */

    CaseSingle cs{ic, icRecent, red, seq};
    // std::vector<RT *> rtVec = cs.apply(graph, forwardVec, backwardVec, graphIndex);


    // Start at 1; we are searching for the nested RT first, then expanding outwards, both before and after. 
    for (int i = 1; i < forwardVec.size(); i++) {
        auto nFor = forwardVec.at(i);
        auto nForConVec = getDiagonals(nFor, graph);

        // Get the first element connected to nFor, ordered by position.
        if (!nForConVec.empty()) {
            auto nBac = nForConVec.front();

            bool isFail = true;

            // A new element to be constructed if there was a hyper-merge
            Element *newEle = nullptr; 
            // The outer left LTR
            Element *oLeftLtr = nullptr;
            // The outer right LTR
            Element *oRightLtr = nullptr;
            // The nested (inner) left LTR
            Element *nLeftLtr = nullptr;
            // The nested (inner) right LTR
            Element *nRightLtr = nullptr;

            // If there was a hyper-merge, outerEle is equal to the unhyper-merged new Element assigned on the heap
            std::tie(newEle, nLeftLtr, nRightLtr) = buildLtrs(nFor, nBac, graph);


            // If the nested LTR RT candidate passes the checkSingle case, expand outward to find an outer RT.
            if (nLeftLtr != nullptr && nRightLtr != nullptr && cs.checkSingle(nLeftLtr, nRightLtr)) {

                oLeftLtr  = findOuterLtr(nFor, graph);
                oRightLtr = findOuterLtr(nBac, graph);

                if (oLeftLtr != nullptr && oRightLtr != nullptr && checkOuterRt(oLeftLtr, nLeftLtr, nRightLtr, oRightLtr)) {
                    // Building strings for Identity check
                    std::string outerInteriorSeq  = seq->substr(oLeftLtr->getEnd(), nLeftLtr->getStart() - oLeftLtr->getEnd());
                    outerInteriorSeq += seq->substr(nRightLtr->getEnd(), oRightLtr->getStart() - nRightLtr->getEnd());
                    std::string innerInteriorSeq  = seq->substr(nLeftLtr->getEnd(), nRightLtr->getStart() - nLeftLtr->getEnd());

                    if (areRecentSeqSame(outerInteriorSeq, innerInteriorSeq)) {
                        isFail = false;

                        nLeftLtr = nLeftLtr == newEle? nLeftLtr : new Element{*nLeftLtr};
                        nRightLtr = nRightLtr == newEle? nRightLtr : new Element{*nRightLtr};
                        rtVec.push_back(new RTComplete{oLeftLtr, oRightLtr, name + "Outer", rank, graphIndex});
                        rtVec.push_back(new RTComplete(nLeftLtr, nRightLtr, name + "Inner", rank, graphIndex));
                    }
                }
            }


            // If the case never converges, delete any element that was created ad hoc
            if (isFail) {
                if (newEle != nullptr) {
                    delete newEle;
                }
                if (oLeftLtr != nullptr) {
                    delete oLeftLtr;
                }
                if (oRightLtr != nullptr) {
                    delete oRightLtr;
                }
            }

        }
    }
}

/**
 * Assuming we are dealing with a recently nested element. 
 * Given one LTR of the nested (inner) element, expand outwards to find the outer (trimmed) LTR.
 * It works for both left and right LTRs.
 * If this is the case, return a new element for the outer LTR.
 */
Element * CaseRecent::findOuterLtr(Element *e, DirectedGraph<Element> &graph) {
    Element *r = nullptr;
    
    // Get the vertical element
    auto v = getVertical(e, graph);

    if (v != nullptr) {
        // get the closest diagonal element
        auto conVec = getDiagonals(v, graph);
        if (!conVec.empty()) {
            // This line is a bit greedy
            Element *c = e->getIsForward() ? conVec.back() : conVec.front();

            if (graph.isBidirectional(*v, *c)) {
                // Overlappnig stretches between e and v
                std::vector<Stretch*> overlapVec = v->getOverlapStretchVec(e);

                int eleStart = -1;
                int eleEnd   = -1;
                std::vector<int> heightVec;

                for (auto s : overlapVec) {
                    // Make a trimmed strech based on the overlap
                    int start = std::max(s->getStart(), e->getStart());
                    int end   = std::min(s->getEnd(), e->getEnd());
                    Stretch cutStretch{start, end, s->getMark(), v->getIsForward()};
                    cutStretch.setMedianHeight(s->getMedianHeight());

                    // Build a hypothetical stretch
                    Stretch hypo = cutStretch.buildMatch();
                    if (c->calcOverlap(hypo) > 0){
                        int hypoStart = hypo.getStart();
                        if (hypoStart < 0) {
                            hypoStart = 0;
                        }
                        int hypoEnd = hypo.getEnd();
                        if (hypoEnd > seq->length()) {
                            hypoEnd = seq->length();
                        }
                        eleStart = eleStart == -1? hypoStart: std::min(hypoStart, eleStart);
                        eleEnd   = std::max(hypo.getEnd(), eleEnd);
                        heightVec.push_back(hypo.getMedianHeight());
                    }
                }
                assert( (eleStart != -1 && eleEnd != -1) || (eleStart == -1 && eleEnd == -1));

                bool order = e->getIsForward() ? e->getStart() > eleEnd : eleStart > e->getEnd();

                if (eleStart != -1 && eleEnd != -1 && order) {
                    r = new Element{eleStart, eleEnd, int(LtrUtility::calcMedian(heightVec)), e->getIsForward() };
                }
            }
        }
    }

    return r;
}

/**
 * Generate a hypothetical RT from a suppossed recently nested element and check if it is a valid single RT.
 */
bool CaseRecent::checkOuterRt(Element *leftLtr, Element *nestedLeftLtr, Element *nestedRightLtr, Element *rightLtr) {
    assert (isAfter(nestedLeftLtr, leftLtr));
    assert (isAfter(nestedRightLtr, nestedLeftLtr));
    assert (isAfter(rightLtr, nestedLeftLtr));


    CaseSingle cs{ic, icRecent, red, seq};
    
    bool r = false;
    
    bool isGap = leftLtr->calcGap(*nestedLeftLtr) > 0 || nestedRightLtr->calcGap(*rightLtr) > 0;

    if (isGap) {
        auto s1 = 0;
        auto s2 = s1 + leftLtr->getSize();
        auto s3 = s2 + nestedLeftLtr->getStart() - leftLtr->getEnd();
        auto s4 = s3 + nestedRightLtr->getEnd() - nestedLeftLtr->getStart();
        auto s5 = s4 + rightLtr->getStart() - nestedRightLtr->getEnd();

        std::vector<int> rtScoreVec = scoreSeq(leftLtr->getStart(), rightLtr->getEnd());

        std::vector<int> leftLtrScoreVec{rtScoreVec.begin() + s1, rtScoreVec.begin() + s2};

        std::vector<int> interiorScoreVec{rtScoreVec.begin() + s2, rtScoreVec.begin() + s3};
        interiorScoreVec.insert(interiorScoreVec.end(), rtScoreVec.begin() + s4, rtScoreVec.begin() + s5);
        
        std::vector<int> rightLtrScoreVec{rtScoreVec.begin() + s5, rtScoreVec.end()};

        r = cs.checkSingle(leftLtrScoreVec, interiorScoreVec, rightLtrScoreVec);
    }

    return r;
}