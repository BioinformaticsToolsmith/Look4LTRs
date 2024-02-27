/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Element.cpp
 *
 *  Created on: Oct 13, 2022
 *      Author: Anthony B. Garza.
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

#include "Element.h"

Element::Element(Stretch &s)
{
    this->start = s.getStart();
    this->end = s.getEnd();
    this->isForward = s.getIsForward();
    pushStretch(s);
}

Element::Element(const Element &other) {
    this->start = other.start;
    this->end = other.end;
    this->isForward = other.isForward;
    for (auto ptr : other.getStretchVec()) {
        stretchVec.push_back(ptr);
    }
}


Element::Element(Element &elementOne, Element &elementTwo)
{
    this->start = std::min(elementOne.getStart(), elementTwo.getStart());
    this->end = std::max(elementOne.getEnd(), elementTwo.getEnd());
    this->isForward = elementOne.getIsForward();
    for (auto ptr : elementOne.getStretchVec())
    {
        stretchVec.push_back(ptr);
    }
    for (auto ptr : elementTwo.getStretchVec())
    {
        stretchVec.push_back(ptr);
    }
}

Element::Element(int start, int end, int height, bool isForward) {
    this->start = start;
    this->end = end;
    this->isForward = isForward;
    Stretch *s = new Stretch{start, end, Stretch::K, isForward};
    s->setMedianHeight(height);
    pushStretch(*s);
    heapVec.push_back(s);
}

Element::Element(std::vector<Stretch*> sVec) {
    assert (!sVec.empty());
    std::sort(sVec.begin(), sVec.end(), [](Stretch *s1, Stretch *s2) {
        return s1->getStart() < s2->getStart();
    });
    this->start = sVec.front()->getStart();
    this->end = sVec.back()->getEnd();
    this->isForward = sVec.front()->getIsForward();
    for (auto s : sVec) {
        assert (s->getIsForward() == isForward);
        pushStretch(*s);
    }
}



Element::~Element() {
    for (auto ptr : heapVec) {
        delete ptr;
    }
}

void Element::merge(Stretch &s)
{
    if (s.getStart() < start)
    {
        start = s.getStart();
    }

    if (s.getEnd() > end)
    {
        end = s.getEnd();
    }

    pushStretch(s);
}

void Element::merge(Element &e) {
    assert (isForward == e.getIsForward());

    start = e.getStart() < start? e.getStart() : start;
    end = e.getEnd() > end? e.getEnd() : end;

    for (auto s : e.getStretchVec()) {
        stretchVec.push_back(s);
    }
}


void Element::pushStretch(Stretch &s)
{
    stretchVec.push_back(&s);
}

// [OK]
bool Element::isOverlap(Element &e)
{
    int overlap = calcOverlap(start, end, e.getStart(), e.getEnd());
    return overlap > 0 ? true : false;
}

// [OK]
int Element::calcOverlap(Stretch &s)
{
    return calcOverlap(start, end, s.getStart(), s.getEnd());
}

// [OK]
int Element::calcOverlap(Element &e)
{
    return calcOverlap(start, end, e.getStart(), e.getEnd());
}

int Element::calcOverlap(const Element &e) {
    return calcOverlap(start, end, e.getStart(), e.getEnd());
}

int Element::calcOverlap(int s, int e)
{
    return calcOverlap(this->start, this->end, s, e);
}


int Element::calcGap(Element &e) {
    int gap;
    int overlap = calcOverlap(e);
    if (overlap == 0) {
        gap = std::max(start, e.getStart()) - std::min(end, e.getEnd());
    }
    else {
        gap = 0;
    }
    return gap;
}

void Element::extendEnd(int amount) {
    end += amount;
}

void Element::extendStart(int amount) {
    start -= amount;
}

int Element::getStart() const
{
    return start;
}

void Element::setStart(int start)
{
    assert(start < end);
    assert(start >= 0);

    this->start = start;
}

int Element::getEnd() const
{
    return end;
}

void Element::setEnd(int end)
{
    assert(end > start);
    assert(end > 0);

    this->end = end;
}

int Element::getSize() const
{
    return end - start;
}

bool Element::getIsForward() const
{
    return isForward;
}

int Element::getMedianHeight() const {
    // Key is the height, size is the value
    int r = -1; // value to return

    std::vector<Stretch *> sVec{stretchVec.begin(), stretchVec.end()};
    std::sort(sVec.begin(), sVec.end(), [](Stretch *a, Stretch *b) {return a->getMedianHeight() < b->getMedianHeight();});

    int total = 0;
    std::vector<int> cumSum (sVec.size(), 0);
    for (int i = 0; i < sVec.size(); i++) {
        total += sVec[i]->getSize();
        cumSum[i] = total;
    }

    if (total % 2 == 1) {
        int med = (total + 1) / 2;
        for (int i = 0; i < cumSum.size(); i++) {
            if (cumSum[i] >= med) {
                r = sVec[i]->getMedianHeight();
                break;
            }
        }
    }
    else {
        int med = total / 2;
        for (int i = 0; i < cumSum.size(); i++) {
            if (cumSum[i] >= med) {
                r = cumSum[i] == med? (sVec[i]->getMedianHeight() + sVec[i + 1]->getMedianHeight()) / 2 : sVec[i]->getMedianHeight();
                break;
            }
        }
    }

    assert(r != -1);

    return r;

}

std::vector<Stretch *> Element::getOverlapStretchVec(Element *e) const {
    std::vector<Stretch *> r;
    for (auto s : stretchVec) {
        if (e->calcOverlap(*s) > 0) {
            r.push_back(s);
        }
    }
    // else {
    //     for (auto s : stretchVec) {
    //         Stretch hypoStretch = s->buildMatch();
    //         Stretch *match = new Stretch{hypoStretch};
    //         if (e->calcOverlap(*match) > 0) {
    //             r.push_back(match);
    //         }
    //         else {
    //             delete match;
    //         }
    //     }
    // }


    return r;
}


Stretch Element::buildStretch() {
    Stretch r{start, end, Stretch::K, isForward};
    r.setMedianHeight(getMedianHeight());
    return r;
}

// std::pair<int, int> Element::findMissingMatch(Element *e) {
//     int minStart = -1;
//     int maxEnd = -1;

//     // Getting the overlapping stretches
//     auto overlapVec = getOverlapStretchVec(e, true);
//     std::unordered_set<int> overlapSet;
//     for (auto s : overlapVec) {
//         overlapSet.insert(s->getStart());
//     }

//     // Getting the stretches that are not overlapping and finding the new start and ends
//     for (auto s : stretchVec) {
//         Stretch match = s->buildMatch();
//         if (overlapSet.count(match.getStart()) == 1) {
//             continue;
//         }

//         if (match.getEnd() < match.getStart()) {
//             minStart = minStart == -1? match.getStart() : std::min(minStart, match.getStart());
//         }
//         else if (match.getStart() > e->getEnd()){
//             maxEnd = std::max(maxEnd, match.getEnd());
//         }
//     }

//     for (auto s : overlapVec) {
//         delete s;
//     }
//     std::pair<int, int> r{minStart, maxEnd};
    
//     return r;
// }

std::vector<Stretch *> Element::getStretchVec() const
{
    std::vector<Stretch *> r{stretchVec.begin(), stretchVec.end()};
    std::sort(r.begin(), r.end(), [](Stretch *s1, Stretch *s2)
            { return s1->getStart() < s2->getStart(); }); 
    return r;
}



bool operator<(const Element &e1, const Element &e2)
{
    return (e1.getStart() < e2.getStart());
}

std::ostream &operator<<(std::ostream &os, const Element &n)
{
    os << n.getStart() << ":" << n.getEnd();
    for (auto x : n.getStretchVec())
    {
        os << " [" << x->getMedianHeight() << "," << x->getStart() << "," << x->getEnd() << "]";
    }
    if (n.getIsForward() == true)
    {
        os << " +";
    }
    else
    {
        os << " -";
    }

    return os;
}

std::pair<int, int> Element::getHypoStarts(Element *e) const {
    int s1 = std::numeric_limits<int>::max();
    int s2 = s1;
    for (auto s : stretchVec) {
        Stretch hypoStretch = s->buildMatch();
        if (e->calcOverlap(hypoStretch) > 0) {
            auto overlap1 = getOverlap(*s);
            auto overlap2 = e->getOverlap(hypoStretch);
            if (overlap1.first < s1){
                s1 = overlap1.first;
                s2 = overlap2.first;                
            }            
        }
    }
    assert(s1 >= 0);
    assert(s2 >= 0);
    // assert(s1 < s2);

    return std::make_pair(s1, s2);
}

std::pair<int, int> Element::getOverlap(Element &e) {
    return getOverlap(start, end, e.getStart(), e.getEnd());
}

std::pair<int, int> Element::getOverlap(Stretch &s) {
    return getOverlap(start, end, s.getStart(), s.getEnd());
}
