/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * CaseMatcher
 *
 *  Created on: Oct 27, 2022
 *      Author: Anthony B. Garza.
 * Reviewer:
 *   Purpose: Abstract class; extended classes are applied to element matching and make decisions based
 *              on their own individual criteria
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

#include "DirectedGraph.h"
#include "Element.h"
#include "RT.h"
#include "../red/Red.h"
#include "../IdentityCalculator.h"
#include "LtrParameters.h"
#include "LtrUtility.h"


#include "LocalAlignment.h"


#include <vector>
#include <algorithm>
#include "assert.h"
#include <unordered_map>
#include <string>

class CaseMatcher
{
public:

    /**
     * Constructor
     */
    CaseMatcher(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq);

    /**
     * Methods
     */
    virtual void apply(DirectedGraph<Element> &graph, std::vector<Element *> &forwardVec, std::vector<Element *> &backwardVec, int graphIndex) = 0;

    bool isReciprocal(Element *left, Element *right, DirectedGraph<Element>& graph);

    Element* isStretchReciprocal(Element *left, Element *right, DirectedGraph<Element>& graph, bool &delLeft, bool &delRight);

    /**
     * Getters
     */
    std::vector<RT *> getRTVec() const;

    // bool hasRTs() const;

    int getRank() const;
    std::string getName() const;




protected:
    /**
     * Variables
     */
    std::vector<RT *> rtVec;
    IdentityCalculator<int32_t> &ic;
    IdentityCalculator<int32_t> &icRecent;
    Red &red;
    const std::string *seq;

    // Rank is used to determine which analysis case to use
    int rank;
    std::string name;
    /**
     * Methods
     */

    // Is the forward element the same as the backwards element? Check if bidirectional with weight 0
    bool isSame(Element *elementOne, Element *elementTwo, DirectedGraph<Element> &graph);

    // Does the element have overlap?
    // bool hasOverlap(Element *e, DirectedGraph<Element> &graph);

    // bool isOverlap(std::vector<Element *> &eVec, DirectedGraph<Element> &graph);

    bool checkLength(std::vector<Element *> &eVec);

    bool checkLength(Element *e);

    // This method calculates the length of a complete RT elmeent and adjusts the end by k - 1
    int adjustSize(int s, int e);

    // Is e1 after e2 location wise? i.e., e1 is to the right of e2 AND there is a gap between.
    bool isAfter(Element *e1, Element *e2);

    // If the element has a vertical connect, get it; otherwise, return a nullptr
    Element *getVertical(Element *e, DirectedGraph<Element> &graph);

    // Build an element with a vertical(overlapping) element if exists
    Element * buildElement(Element *e, DirectedGraph<Element> &graph, bool &isDelete);

    std::vector<Element*> getDiagonals(Element *e, DirectedGraph<Element> & graph);

    std::vector<Element*> getHorizontals(Element *e, DirectedGraph<Element> &graph, bool forward);

    std::vector<Element*> getHorizontalsAcross(Element *e, DirectedGraph<Element> &graph, bool forward);

    std::tuple<Element*, Element*, Element*> buildLtrs(Element *f, Element *b, DirectedGraph<Element> &graph);

    std::tuple<Element*, Element*> assignLtrs(Element *built, Element *left, Element *right);

    std::vector<int> scoreSeq(int s, int e);

    std::string substrInterior(Element * e1, Element *e2, std::string *seq);

    bool areSeqSame(std::string &seq1, std::string &seq2);

    bool areRecentSeqSame(std::string &seq1, std::string &seq2);

    Element * retrieveFirst(Element *e1, Element *e2);

    double calcSizeRatio(int s1, int e1, int s2, int s3);

    double calcSizeRatio(int size1, int size2);
};