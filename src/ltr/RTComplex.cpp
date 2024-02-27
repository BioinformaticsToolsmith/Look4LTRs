/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * RTComplex
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

#include "RTComplex.h"

RTComplex::RTComplex(std::vector<Element*> _ltrVec, std::string caseType, int caseRank, int graphGroup) : ltrVec{_ltrVec}
{
    this->ltrVec = ltrVec;
    this->caseType = caseType;
    this->caseRank = caseRank;
    this->graphGroup = graphGroup;
}

RTComplex::~RTComplex() {
    for (auto ltr : getLTRVec()) {
        delete ltr;
    }
}

int RTComplex::getStart() const
{
    return ltrVec.front()->getStart();
}

int RTComplex::getEnd() const
{
    return ltrVec.back()->getEnd();
}

int RTComplex::getSize() const {
    return ltrVec.back()->getEnd() - ltrVec.front()->getStart();
}

int RTComplex::getIntSize(bool withNested = true) const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular interior." << std::endl;
    throw std::exception();
}

std::string RTComplex::getIntSeq(std::string &seq, bool withNested = true) const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular interior." << std::endl;
    throw std::exception();
}


Range RTComplex::getRange() const {
    Range r;
    for (auto ltr : ltrVec) {
        r.push_back(std::pair<int, int>{ltr->getStart(), ltr->getEnd()});
    }

    return r;
}

bool RTComplex::getIsRC() const {
    std::cerr << "Unsupported operation: A complex LTR RT can be a mix of RCs and non-RCs." << std::endl;
    throw std::exception();
}

bool RTComplex::getIsTSDExist() const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular TSD." << std::endl;
    throw std::exception();
}

const Element* RTComplex::getLeftLTR() const
{
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular Left LTR." << std::endl;
    throw std::exception();
    }

const Element* RTComplex::getRightLTR() const
{
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular right LTR." << std::endl;
    throw std::exception();
}

std::vector<Element*> RTComplex::getLTRVec() const {
    return ltrVec;
}

const std::set<RT*> RTComplex::getNestSet() const {
    std::cerr << "Unsupported operation: A complex LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

const std::set<RT*> RTComplex::getOuterSet() const {
    std::cerr << "Unsupported operation: A complex LTR RT can not be nested." << std::endl;
    throw std::exception();}

int RTComplex::getPPTStart() const {
    std::cerr << "Unsupported operation: A complex LTR RT has multiple possible PPTs." << std::endl;
    throw std::exception();
}

int RTComplex::getPPTEnd() const {
    std::cerr << "Unsupported operation: A complex LTR RT has multiple possible PPTs." << std::endl;
    throw std::exception();
}

std::pair<int, int> RTComplex::getLeftTSD() const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular TSD." << std::endl;
    throw std::exception();
}

std::pair<int, int> RTComplex::getRightTSD() const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular TSD." << std::endl;
    throw std::exception();
}

double RTComplex::getIdentityScore() const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have an identity score." << std::endl;
    throw std::exception();
}

void RTComplex::setIsRC(bool isRC) {
    std::cerr << "Unsupported operation: Unable to set an complex LTR RT to be RC." << std::endl;
    throw std::exception();
}

void RTComplex::setPPT(int start, int end) {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular PPT." << std::endl;
    throw std::exception();
}

void RTComplex::setTSD(int leftStart, int leftEnd, int rightStart, int rightEnd) {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular TSD." << std::endl;
    throw std::exception();
}

void RTComplex::setIdentityScore(double identityScore) {
    std::cerr << "Unsupported operation: Unable to set an complex LTR RT to have an identity score." << std::endl;
    throw std::exception();
}


bool RTComplex::hasRightLTR() const
{
    return false;
}

bool RTComplex::hasLeftLTR() const {
    return false;
}

void RTComplex::nest(RT* rt) {
    std::cerr << "Unsupported operation: A complex LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

void RTComplex::removeNest(RT* rt) {
    std::cerr << "Unsupported operation: A complex LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

bool RTComplex::hasNest() const {
    std::cerr << "Unsupported operation: A complex LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

bool RTComplex::isNested() const {
    std::cerr << "Unsupported operation: A complex LTR RT can not be nested." << std::endl;
    throw std::exception();
}

bool RTComplex::couldNest(RT *rt) const {
    return false;
}

void RTComplex::addOuter(RT *rt) {
    std::cerr << "Unsopported operation: A complex LTR RT can not be nested." << std::endl;
    throw std::exception();
}

void RTComplex::removeOuter(RT *rt) {
    std::cerr << "Unsopported operation: A complex LTR RT can not be nested." << std::endl;
    throw std::exception();
}


void RTComplex::extend(int k, bool isForward) {
    for (auto ltrPtr : ltrVec) {
        if (isForward) {
            ltrPtr->extendEnd(k);
        }
        else {
            ltrPtr->extendStart(k);
        }
    }
}

void RTComplex::expand(int pos, int length) {

    for (auto ltrPtr : ltrVec) {
        if (pos <= ltrPtr->getEnd()) {
            ltrPtr->setEnd(ltrPtr->getEnd() + length);
        }
        if (pos < ltrPtr->getStart()) {
            ltrPtr->setStart(ltrPtr->getStart() + length);
        }
    }

}

void RTComplex::push(int p) {

    for (auto ltrPtr : ltrVec) {
        ltrPtr->setEnd(ltrPtr->getEnd() + p);
        ltrPtr->setStart(ltrPtr->getStart() + p);
    }
}


std::ostream &operator<<(std::ostream &os, const RTComplex &r)
{
    auto rangeVec = r.getRange();
    for (int i = 0; i < rangeVec.size() - 1; i++) {
        os << rangeVec.at(i).first << " " << rangeVec.at(i).second << " ";
    }
    os << rangeVec.back().first << " " << rangeVec.back().second;
    return os;
}