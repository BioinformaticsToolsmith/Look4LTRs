/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Element.h
 *
 *  Created on: Oct 13, 2022
 *      Author: Anthony B. Garza.
 *    Reviewer: Hani Z. Girgis.
 *
 * Purpose: Contains the start and ends of an element in a chromosome
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

#include "Stretch.h"

#include <assert.h>
#include <vector>
#include <memory>
#include <algorithm>
#include <unordered_set>
#include <limits>

class Element
{

public:
    /**
     * Constructors
     */

    Element (const Element &other);

    // This constructor takes a stretch's start and end as well as adds the median height to its heightVec
    Element(Stretch &s);

    // This constructor takes a vector of stretches
    Element(std::vector<Stretch*> sVec);

    // This constructor merges two elements to build a new element
    Element(Element &elementOne, Element &elementTwo);

    // This constructor builds an element from a start, end, height, and direction
    Element(int start, int end, int height, bool isForward);

    ~Element();

    /**
     * Getter and Setters
     */
    int getStart() const;
    void setStart(int start);

    int getEnd() const;
    void setEnd(int end);

    int getSize() const;

    bool getIsForward() const;

    int getMedianHeight() const;

    std::vector<Stretch *> getOverlapStretchVec(Element *e) const;

    std::pair<int, int> getHypoStarts(Element *e) const;

    std::vector<Stretch *> getStretchVec() const;

    /**
     * Methods
     */

    // Update the start and ends of the element as well as the stretchVec according to the given stretch
    void merge(Stretch &s);

    // Update the start and ends of the element according to another element
    void merge(Element &e);

    // Push the stretch into the stretchVec
    void pushStretch(Stretch &s);

    // Does the element overlap with the given element?
    bool isOverlap(Element &e);

    // Calculate the amount of overlap between the element and a given stretch
    int calcOverlap(Stretch &s);

    // Calculate the amount of overlap between the element and a given element
    int calcOverlap(Element &e);
    int calcOverlap(const Element &e);

    std::pair<int, int> getOverlap(Element &e);
    std::pair<int, int> getOverlap(Stretch &s);

    // Calculate the amount of overlap between the element and a given start and end
    int calcOverlap(int s, int e);

    // Calculate the gap between two elements; if overlap, provide the overlap as a negative number
    int calcGap(Element &e);

    void extendEnd(int amount);
    void extendStart(int amount);

    Stretch buildStretch();

    std::pair<int, int> findMissingMatch(Element *e);

    // Less than operator for comparison and sorting

private:
    /**
     * Variables
     */
    // start is inclusive, starting at 0
    int start;
    // end is exclusive
    int end;

    // is the element pointing forward?
    bool isForward;

    // Vector of pointers to stretches
    std::vector<Stretch *> stretchVec;

    std::vector<Stretch *> heapVec;

    inline int calcOverlap(int s1, int e1, int s2, int e2)
    {
        int overlap = std::min(e1, e2) - std::max(s1, s2);
        return overlap > 0 ? overlap : 0;
    }

    std::pair<int, int> getOverlap(int s1, int e1, int s2, int e2) {
        return std::make_pair(std::max(s1, s2), std::min(e1, e2));        
    }
};

bool operator<(const Element &e1, const Element &e2);

std::ostream &operator<<(std::ostream &os, const Element &n);