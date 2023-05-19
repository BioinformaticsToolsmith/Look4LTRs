/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * CaseSolo
 *
 *  Created on: Oct 27, 2022
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

#include "CaseSolo.h"

CaseSolo::CaseSolo(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq) : CaseMatcher(_ic, _icRecent, _red, _seq)
{
    name = "Solo";
    rank = 150;
}

void CaseSolo::apply(DirectedGraph<Element> &graph, std::vector<Element *> &forwardVec, std::vector<Element *> &backwardVec, int graphIndex)
{

    /**
     * Simple Single RT with a matching solo LTR on either side (but not both), three solo LTRS or a nested solo RT
     *
     * >  >-->   or   >-->  >  or   >-->-->
     *
     * O  O
     *  \ |\
     *   \| \
     *    O  O
     *
     * Either the top-left node, middle node, or the bottom-right node may represent the solo LTR.  Additionally, the top-left node may
     * match all bottom nodes and the bottom-right node may match all top nodes.
     *
     * Check for two forward and two backward nodes.  Ensure that the first forward has no overlap with the first backward postion-wise
     * (an LTR shouldn't overlap with another LTR) as well as for the second forward and second backward.
     * Then, check if the second forward is overlapped with the first forward (They should be the same LTR).
     *
     * To determine which side the solo LTR is on, Red scores should be used on the areas between the nodes. If between two nodes
     * there is a high Red score, then that area is most likely an interior.  If a low Red score, then one of the nodes must be
     * the solo LTR. If both the interiors together have a high red score, then the solo LTR is between.
     */

    CaseSingle cs{ic, icRecent, red, seq};

    try {
    // minimum required for this structure
    if (forwardVec.size() > 1 && backwardVec.size() > 1)
    {
        bool isComplete = false;
        for (auto fElePtr : forwardVec)
        {
            if (isComplete)
            {
                break;
            }

            // Look Forward
            for (auto bElePtr : getDiagonals(fElePtr, graph))
            {
                if (isComplete)
                {
                    break;
                }

                // building left and middle LTRs based on reciprocal coverage and/or hyper-merging
                if (graph.isBidirectional(*fElePtr, *bElePtr) && isReciprocal(fElePtr, bElePtr, graph))
                {
                    auto vert = getVertical(bElePtr, graph);
                    if (vert != nullptr)
                    {
                        for (auto finalElePtr : getDiagonals(vert, graph))
                        {
                            if (graph.isBidirectional(*vert, *finalElePtr) && isReciprocal(vert, finalElePtr, graph))
                            {
                                Element *middleLtr = new Element{*bElePtr, *vert};

                                if (fElePtr->getEnd() <= middleLtr->getStart() && middleLtr->getEnd() <= finalElePtr->getStart())
                                {
                                    bool isLeft = cs.checkSingle(fElePtr, middleLtr);
                                    bool isRight = cs.checkSingle(middleLtr, finalElePtr);
                                    bool isAll = checkNestedSolo(fElePtr, middleLtr, finalElePtr);

                                    if (isLeft)
                                    {
                                        rtVec.push_back(new RTComplete{new Element{*fElePtr}, middleLtr, name + "Single", rank, graphIndex});
                                        rtVec.push_back(new RTSolo{new Element{*finalElePtr}, name + "LTR", rank, graphIndex});
                                    }
                                    if (isRight)
                                    {
                                        rtVec.push_back(new RTSolo{new Element{*fElePtr}, name + "LTR", rank, graphIndex});
                                        rtVec.push_back(new RTComplete{middleLtr, new Element{*finalElePtr}, name + "Single", rank, graphIndex});
                                    }
                                    if (isAll)
                                    {
                                        rtVec.push_back(new RTComplete{new Element{*fElePtr}, new Element{*finalElePtr}, name + "Single", rank, graphIndex});
                                        rtVec.push_back(new RTSolo{middleLtr, name + "LTR", rank, graphIndex});
                                    }
                                    if (isLeft || isRight || isAll)
                                    {
                                        isComplete = true;
                                        break;
                                    }
                                }
                                delete middleLtr;
                            }
                        }
                    }
                }
            }
        }
    }
    } catch (...) {
        std::cout << "Solo" << std::endl;
        throw std::exception();
    }
}

bool CaseSolo::checkNestedSolo(Element *leftLtr, std::vector<Element *> nestVec, Element *rightLtr)
{
    bool r = false;
    CaseSingle cs{ic, icRecent, red, seq};
    Element *last = leftLtr;
    bool properOrder = true;
    for (auto nest : nestVec)
    {
        if (nest->getStart() < last->getEnd())
        {
            properOrder = false;
            break;
        }
    }

    if (properOrder)
    {

        std::vector<int> rtScoreVec = scoreSeq(leftLtr->getStart(), rightLtr->getEnd());

        auto s1 = 0;
        auto s2 = s1 + leftLtr->getSize();
        auto s3 = s2 + rightLtr->getStart() - leftLtr->getEnd();
        int start = leftLtr->getStart();

        std::vector<std::vector<int>> ltrScoreVec = {std::vector<int>{rtScoreVec.begin() + s1, rtScoreVec.begin() + s2}};
        std::vector<int> interiorScoreVec;

        last = leftLtr;
        for (auto nest : nestVec)
        {
            ltrScoreVec.push_back(std::vector<int>{rtScoreVec.begin() + nest->getStart() - start,
                                                   rtScoreVec.begin() + nest->getEnd() - start});
            interiorScoreVec.insert(interiorScoreVec.end(), rtScoreVec.begin() + last->getEnd() - start,
                                    rtScoreVec.begin() + nest->getStart() - start);
            last = nest;
        }

        ltrScoreVec.push_back(std::vector<int>{rtScoreVec.begin() + s3, rtScoreVec.end()});
        interiorScoreVec.insert(interiorScoreVec.end(), rtScoreVec.begin() + last->getEnd() - start,
                                rtScoreVec.begin() + rightLtr->getStart() - start);

        bool isGap = interiorScoreVec.size() > 0;
        if (isGap)
        {

            bool failedRatio = false;
            double interiorMedian = LtrUtility::calcMedian(interiorScoreVec);
            for (auto ltrScore : ltrScoreVec)
            {
                double ltrMedian = LtrUtility::calcMedian(ltrScore);
                if (interiorMedian / ltrMedian < LtrParameters::MIN_LTR_INT_RATIO)
                {
                    failedRatio = true;
                    break;
                }
            }

            if (failedRatio == false)
            {
                double interiorPerc = LtrUtility::calcPercent(interiorScoreVec);
                r = interiorPerc >= LtrParameters::MIN_PERC;
            }
        }
    }

    return r;
}

bool CaseSolo::checkNestedSolo(Element *leftLtr, Element *nestedLtr, Element *rightLtr)
{
    bool r = false;
    CaseSingle cs{ic, icRecent, red, seq};

    bool isGap = leftLtr->calcGap(*nestedLtr) > 0 || nestedLtr->calcGap(*rightLtr) > 0;

    if (isGap)
    {

        std::vector<int> rtScoreVec = scoreSeq(leftLtr->getStart(), rightLtr->getEnd());

        auto s1 = 0;
        auto s2 = s1 + leftLtr->getSize();
        auto s3 = s2 + rightLtr->getStart() - leftLtr->getEnd();
        auto nStart = s2 + nestedLtr->getStart() - leftLtr->getEnd();
        auto nEnd = nStart + nestedLtr->getSize();

        std::vector<int> leftLtrScoreVec{rtScoreVec.begin() + s1, rtScoreVec.begin() + s2};
        std::vector<int> interiorScoreVec{rtScoreVec.begin() + s2, rtScoreVec.begin() + nStart};
        interiorScoreVec.insert(interiorScoreVec.end(), rtScoreVec.begin() + nEnd, rtScoreVec.begin() + s3);
        std::vector<int> rightLtrScoreVec{rtScoreVec.begin() + s3, rtScoreVec.end()};
        std::vector<int> nestedLtrScoreVec{rtScoreVec.begin() + nStart, rtScoreVec.begin() + nEnd};

        double interiorMedian = LtrUtility::calcMedian(interiorScoreVec);
        double nestedLtrMedian = LtrUtility::calcMedian(nestedLtrScoreVec);
        r = cs.checkSingle(leftLtrScoreVec, interiorScoreVec, rightLtrScoreVec);
    }

    return r;
}