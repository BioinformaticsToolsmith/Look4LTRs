#include "RTComplete.h"

RTComplete::RTComplete(Element *leftLTR, Element *rightLTR, std::string caseType, int caseRank, int graphGroup)
{
    assert(leftLTR != nullptr);
    assert(rightLTR != nullptr);
 
    assert(leftLTR->getStart() < rightLTR->getStart());
    assert(leftLTR->getEnd() <= rightLTR->getStart());

    this->leftLTR = leftLTR;
    this->rightLTR = rightLTR;
    this->caseType = caseType;
    this->caseRank = caseRank;
    this->graphGroup = graphGroup;

    isRC = false;
    PPTStart = -1;
    PPTEnd = -1;
    isTSDExist = false;
    identityScore = -1;
}

RTComplete::~RTComplete() {
    for (auto ptr : nestSet) {
        ptr->removeOuter(this);
    }
    delete leftLTR;
    delete rightLTR;
}

int RTComplete::getStart() const
{
    return leftLTR->getStart();
}

int RTComplete::getEnd() const
{
    return rightLTR->getEnd();
}

int RTComplete::getSize() const {
    return rightLTR->getEnd() - leftLTR->getStart();
}

int RTComplete::getIntSize(bool withNested = true) const {
    int r = rightLTR->getStart() - leftLTR->getEnd();
    // Just get the length of the interior including an nested
    if (!withNested) {
        for (auto &j : nestSet) {
            r -= j->getSize();
        }
    }
    return r;
}

std::string RTComplete::getIntSeq(std::string &seq, bool withNested = true) const {
    std::string r{};

    Range range = getRange();
    for (int i = 1; i < range.size() - 1; i++) {
        r += seq.substr(range[i].first, range[i].second - range[i].first);
    }

    return r;
}


Range RTComplete::getRange() const{
    Range r;
    int start = leftLTR->getStart();
    int end = leftLTR->getEnd();
    r.push_back( std::pair<int, int>{start, end} );

    start = leftLTR->getEnd();
    for (auto &j : nestSet) {
        end = j->getStart();
        if (start < end) {
            r.push_back( std::pair<int, int>{start, end} );
        }
        start = j->getEnd();
    }
    end = rightLTR->getStart();
    r.push_back( std::pair<int, int>{start, end} );

    start = rightLTR->getStart();
    end = rightLTR->getEnd();
    r.push_back( std::pair<int, int>{start, end} );


    return r;
}

bool RTComplete::getIsRC() const {
    return isRC;
}

bool RTComplete::getIsTSDExist() const {
    return isTSDExist;
}

const Element*  RTComplete::getLeftLTR() const
{
    return leftLTR;
}

const Element* RTComplete::getRightLTR() const
{
    return rightLTR;
}

std::vector<Element*> RTComplete::getLTRVec() const {
    std::vector<Element *> r = {leftLTR, rightLTR};
    return r;
}

const std::set<RT*> RTComplete::getNestSet() const {
    return nestSet;
}

const std::set<RT*> RTComplete::getOuterSet() const {
    return outerSet;
}

int RTComplete::getPPTStart() const {
    return PPTStart;
}

int RTComplete::getPPTEnd() const {
    return PPTEnd;
}

std::pair<int, int> RTComplete::getLeftTSD() const {
    return std::make_pair(leftTSDStart, leftTSDEnd);
}

std::pair<int, int> RTComplete::getRightTSD() const {
    return std::make_pair(rightTSDStart, rightTSDEnd);
}

double RTComplete::getIdentityScore() const {
    assert(identityScore >= 0 && identityScore <= 1);
    return identityScore;
}


void RTComplete::setIsRC(bool isRC = true) {
    this->isRC = isRC;
}

void RTComplete::setPPT(int start, int end) {
    assert (start >= 0);
    assert (end > 0);
    assert (start < end);
    assert (start >= leftLTR->getEnd());
    assert (end <= rightLTR->getStart());

    PPTStart = start;
    PPTEnd = end;
}


/**
 * leftStart is the start of the left TSD
 * leftEnd is the end of the left TSD
 * rightStart is the start of the right TSD
 * rightEnd is the end of the right TSD
*/
void RTComplete::setTSD(int leftStart, int leftEnd, int rightStart, int rightEnd) {
    assert(leftStart < leftEnd);
    assert(leftEnd < rightStart);
    assert(rightStart < rightEnd);
    assert(leftEnd <= leftLTR->getStart());
    assert(rightStart >= rightLTR->getEnd());

    isTSDExist = true;

    leftTSDStart = leftStart;
    leftTSDEnd = leftEnd;
    rightTSDStart = rightStart;
    rightTSDEnd = rightEnd;

    leftLTR->setStart(leftTSDEnd);
    rightLTR->setEnd(rightTSDStart);
}

void RTComplete::setIdentityScore(double identityScore) {
    this->identityScore = identityScore;
}


bool RTComplete::hasRightLTR() const
{
    return true;
}

bool RTComplete::hasLeftLTR() const
{
    return true;
}

void RTComplete::nest(RT* rt) {
    nestSet.insert(rt);
    rt->addOuter(this);
}

void RTComplete::removeNest(RT* rt) {
    if (nestSet.count(rt) == 1) {
        nestSet.erase(rt);
    }
}

bool RTComplete::hasNest() const {
    return nestSet.size() > 0 ? true : false;
}

bool RTComplete::isNested() const {
    return outerSet.size() > 0 ? true : false;
}


bool RTComplete::couldNest(RT *rt) const {
    return (rt->getStart() > leftLTR->getEnd() && rt->getEnd() < rightLTR->getStart() &&
            rt->hasLeftLTR()) ? true : false;
}

void RTComplete::addOuter(RT *rt) {
    outerSet.insert(rt);
}

void RTComplete::removeOuter(RT *rt) {
    if (outerSet.count(rt) == 1) {
        outerSet.erase(rt);
    }
}

void RTComplete::extend(int k, bool isForward) {
    int lend = leftLTR->getEnd();
    int rstart = rightLTR->getStart();
    if (isForward && lend + k < rstart) {
        leftLTR->extendEnd(k);
        rightLTR->extendEnd(k);
    }
    else if (!isForward && lend < rstart - k) {
        leftLTR->extendStart(k);
        rightLTR->extendStart(k);
    }


}

void RTComplete::expand(int pos, int length) {
    if (pos <= leftLTR->getEnd()) {
        leftLTR->setEnd(leftLTR->getEnd() + length);
    }
    if (pos < leftLTR->getStart()) {
        leftLTR->setStart(leftLTR->getStart() + length);
    }
    if (pos <= rightLTR->getEnd()) {
        rightLTR->setEnd(rightLTR->getEnd() + length);
    }
    if (pos < rightLTR->getStart()) {
        rightLTR->setStart(rightLTR->getStart() + length);
    }

}

void RTComplete::push(int p) {
    int lstart = leftLTR->getStart();
    int lend = leftLTR->getEnd();
    int rstart = rightLTR->getStart();
    int rend = rightLTR->getEnd();

    rightLTR->setEnd(rend + p);
    rightLTR->setStart(rstart + p);
    leftLTR->setEnd(lend + p);
    leftLTR->setStart(lstart + p);
}

std::ostream &operator<<(std::ostream &os, const RTComplete &r)
{
    os << r.getLeftLTR()->getStart() << " " << r.getLeftLTR()->getEnd() << " " << r.getRightLTR()->getStart() << " " << r.getRightLTR()->getEnd();
    return os;
}