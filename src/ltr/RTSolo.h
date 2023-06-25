/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * RTSolo
 *
 *  Created on: X X, 20XX
 *      Author: Anthony B. Garza.
 *      Author: Hani Z. Girgis
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
#include "Element.h"

#include <string>
#include <set>

class RTSolo : public RT
{
private:
    Element *ltr;
    std::set<RT*> outerSet;

public:
    RTSolo(Element *ltr, std::string caseType, int caseRank, int graphGroup);
    ~RTSolo();

    virtual int getStart() const;
    virtual int getEnd() const;
    virtual int getSize() const;
    virtual int getIntSize(bool withNested = true) const;
    virtual std::string getIntSeq(std::string &seq, bool withNested = true) const;
    Range getRange() const;
    virtual bool getIsRC() const;
    virtual bool getIsTSDExist() const;
    virtual const Element *getLeftLTR() const;
    virtual const Element *getRightLTR() const;
    virtual std::vector<Element*> getLTRVec() const override;
    virtual const std::set<RT*> getNestSet() const;
    virtual const std::set<RT*> getOuterSet() const;
    virtual int getPPTStart() const;
    virtual int getPPTEnd() const;
    virtual std::pair<int, int> getLeftTSD() const;
    virtual std::pair<int, int> getRightTSD() const;
    virtual double getIdentityScore() const;


    virtual void setIsRC(bool isRC = true);
    virtual void setPPT(int start, int end);
    virtual void setTSD(int leftStart, int leftEnd, int rightStart, int rightEnd);
    virtual void setIdentityScore(double identityScore);

    virtual bool hasRightLTR() const;
    virtual bool hasLeftLTR() const;

    virtual void nest(RT* rt);
    virtual void removeNest(RT* rt);
    virtual bool hasNest() const;
    virtual bool isNested() const;

    bool couldNest(RT *rt) const;
    void addOuter(RT *rt);
    void removeOuter(RT *rt);

    void extend(int k, bool isForward = true);
    void expand(int pos, int length);
    void push(int p);

};

std::ostream &operator<<(std::ostream &os, const RTSolo &r);
