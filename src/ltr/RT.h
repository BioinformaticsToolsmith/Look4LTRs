/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * RT.h
 *
 *  Created on: Oct 19, 2022
 *      Author: Anthony B. Garza.
 *      Author: Hani Z. Girgis.
 *     Purpose: An interface for complete LTR RTs and solo LTR RTs
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

#include "Element.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include <string>


typedef std::vector<std::pair<int, int>> Range;

class RT
{
protected:
    // If this RT is nested, outerNest points to the RT in which this one is nested in
    RT* outerNest = nullptr; 

    std::string caseType;
    int caseRank;
    int graphGroup;

public:

    /**
     * Getters
     */

    // Get start of the RT
    virtual int getStart() const = 0;
    // Get end of the RT   
    virtual int getEnd() const = 0;   
    // Get size of the RT  
    virtual int getSize() const = 0; 
    // Get size of the interior   
    virtual int getIntSize(bool withNested = true) const = 0; 
    // Get the sequence of the interior
    virtual std::string getIntSeq(std::string &seq, bool withNested = true) const = 0;
    // Get the ranges of starts and ends; 5' start to 5' end, 5'end to nested 1's start, nested 1's end to nested 2's start... nested x's end to 3' start, 3' start to 3' end
    virtual Range getRange() const = 0;
    // Get if the RT is reverse complemented
    virtual bool getIsRC() const = 0;   
    // Get if the RT has a TSD
    virtual bool getIsTSDExist() const = 0;
    // Get the left LTR
    virtual const Element *getLeftLTR() const = 0;
    // Get the right LTR
    virtual const Element *getRightLTR() const = 0;
    // Get all of the LTRs in the RT
    virtual std::vector<Element*> getLTRVec() const = 0;
    // Get the set of nested RTs
    virtual const std::set<RT*> getNestSet() const = 0;
    // Get PPT start
    virtual int getPPTStart() const = 0;
    // Get PPT end
    virtual int getPPTEnd() const = 0;
    // Get TSD start
    virtual std::pair<int, int> getLeftTSD() const = 0;
    // Get TSD end
    virtual std::pair<int, int> getRightTSD() const = 0;
    // Get the Identity score between the two LTRs for only RT COMPLETE
    virtual double getIdentityScore() const = 0;

    // Get the outermost RT that this RT is nested in
    RT* getOuterNest() {
        return outerNest;
    }
    // Get the case type
    std::string getCaseType() {
        return caseType;
    }
    // Get the case rank
    int getCaseRank() {
        return caseRank;
    }
    // Get the graph group (the group of RTs that match each other by LTRs)
    int getGraphGroup() {
        return graphGroup;
    }

    /**
     * Setters
    */
    // Set RC
    virtual void setIsRC(bool isRC = true) = 0;
    // Set the PPT. This is in the interior of the RT
    virtual void setPPT(int start, int end) = 0;
    // Set the TSD.
    virtual void setTSD(int leftStart, int leftEnd, int rightStart, int rightEnd) = 0;
    // Set the Identity score between the two LTRs of a complete RT
    virtual void setIdentityScore(double identityScore) = 0;
    // Rename the Case Type
    void setCaseType(std::string caseType) {
        this->caseType = caseType;
    }
    // Set the outermost RT that this RT is nested in
    void setOuterNest(RT* rt) {
        this->outerNest = rt;
    }
    // Sets the graph group
    void setGraphGroup(int graphGroup) {
        this->graphGroup = graphGroup;
    }

    // Does the LTR RT have a right LTR?
    virtual bool hasRightLTR() const = 0;
    // Does the LTR RT have a left LTR?
    virtual bool hasLeftLTR() const = 0;
    // Does the LTR RT have a TSD?

    /**
     * Nest methods
    */
    // Note that the given RT is nested in this RT
    virtual void nest(RT* rt) = 0;
    // Note that the given RT is no longer nested in this RT
    virtual void removeNest(RT* rt) = 0;
    // Get if this RT has any nested RTs
    virtual bool hasNest() const = 0;
    // Get if the given RT could be nested in this RT
    virtual bool couldNest(RT *rt) const = 0;
    // Get if this RT is nested in another RT
    bool hasOuter() {
        return outerNest != nullptr;
    }


    /**
     * Auxiliary methods
     */
    // extends the LTR RT by k amount
    virtual void extend(int k, bool isForward = true) = 0;
    // Expand the RT outwards from a certain position by length amount
    virtual void expand(int pos, int length) = 0;
    // Push the RT forward a certain amount
    virtual void push(int ps) = 0;
    // Get the overlap between this RT and the given RT
    int calcOverlap(RT *r1) {
        return std::min(getEnd(), r1->getEnd()) - std::max(getStart(), r1->getStart());
    }

};


inline bool operator<(const RT &r1, const RT &r2)
{
    return r1.getStart() < r2.getStart();
}

inline std::ostream &operator<<(std::ostream &os, const RT * r) {
    os << r->getStart() << " " << r->getLeftLTR()->getEnd() << " ";
    if (r->hasRightLTR()) {
        os << r->getRightLTR()->getStart() << " " << r->getEnd();
    }
    return os;
}