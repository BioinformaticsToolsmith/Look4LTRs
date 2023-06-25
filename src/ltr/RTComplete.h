/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * RTComplete
 *
 *  Created on: X X, 20XX
 *      Author: Hani Z. Girgis
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

#include "RT.h"
#include "assert.h"

#include <set>
#include <unordered_set>
#include <algorithm>

typedef std::vector<std::pair<int, int>> Range;


class RTComplete : public RT
{
private:
    Element *leftLTR;
    Element *rightLTR;
    
    bool isRC;
    int PPTStart;
    int PPTEnd;

    bool isTSDExist;
    int leftTSDStart;
    int leftTSDEnd;
    int rightTSDStart;
    int rightTSDEnd;
    double identityScore;

    // std::set<RT*, bool(*)(RT*, RT*)> nestSet{[](RT* a, RT* b){
    //     if (a->getStart() == b->getStart()) {
    //         return a->getEnd() < b->getEnd();
    //     } else {
    //         return a->getStart() < b->getStart();
    //     }
    // }};
    std::set<RT*> nestSet;
    std::set<RT*> outerSet;
public:
    RTComplete(Element *leftLTR, Element *rightLTR, std::string caseType, int caseRank, int graphGroup);
    ~RTComplete();

    int getStart() const;
    int getEnd() const;
    int getSize() const;
    int getIntSize(bool withNested = true) const;
    std::string getIntSeq(std::string &seq, bool withNested = true) const;
    Range getRange() const;
    bool getIsRC() const;
    bool getIsTSDExist() const;
    const Element *getLeftLTR() const;
    const Element *getRightLTR() const;
    std::vector<Element*> getLTRVec() const override;
    const std::set<RT*> getNestSet() const;
    const std::set<RT*> getOuterSet() const;
    int getPPTStart() const;
    int getPPTEnd() const;
    std::pair<int, int> getLeftTSD() const;
    std::pair<int, int> getRightTSD() const;
    double getIdentityScore() const;

    void setIsRC(bool isRC = true);
    void setPPT(int start, int end);
    void setTSD(int leftStart, int leftEnd, int rightStart, int rightEnd);
    void setIdentityScore(double identityScore);


    bool hasRightLTR() const;
    bool hasLeftLTR() const;

    void nest(RT* rt);
    void removeNest(RT* rt);
    bool hasNest() const;
    bool isNested() const;
    bool couldNest(RT *rt) const;
    void addOuter(RT *rt);
    void removeOuter(RT *rt);

    void extend(int k, bool isForward = true);
    void expand(int pos, int length);
    void push(int ps);

};

std::ostream &operator<<(std::ostream &os, const RTComplete &r);
