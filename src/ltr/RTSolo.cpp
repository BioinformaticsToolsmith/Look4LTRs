#include "RTSolo.h"

RTSolo::RTSolo(Element *ltr, std::string caseType, int caseRank, int graphGroup)
{
    assert(ltr != nullptr);

    this->ltr = ltr;
    this->caseType = caseType;
    this->caseRank = caseRank;
    this->graphGroup = graphGroup;
}

RTSolo::~RTSolo() {
    delete ltr;
}

int RTSolo::getStart() const
{
    return ltr->getStart();
}

int RTSolo::getEnd() const
{
    return ltr->getEnd();
}

int RTSolo::getSize() const {
    return ltr->getSize();
}

int RTSolo::getIntSize(bool withNested = true) const {
    std::cerr << "Unsupported operation: A solo LTR RT does not have an interior." << std::endl;
    throw std::exception();
}

std::string RTSolo::getIntSeq(std::string &seq, bool withNested = true) const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have an interior." << std::endl;
    throw std::exception();
}


Range RTSolo::getRange() const {
    Range r{ std::pair<int, int>{ltr->getStart(), ltr->getEnd()} };
    return r;
}

bool RTSolo::getIsRC() const {
    std::cerr << "Unsupported operation: Unable to determine if a solo LTR is RC." << std::endl;
    throw std::exception();
}

bool RTSolo::getIsTSDExist() const {
    std::cerr << "Unsupported operation: A complex LTR RT does not have a singular TSD." << std::endl;
    throw std::exception();
}

const Element* RTSolo::getLeftLTR() const
{
    return ltr;
}

const Element* RTSolo::getRightLTR() const
{
    std::cerr << "Unsupported operation: A solo LTR RT does not have a right LTR." << std::endl;
    throw std::exception();
}

std::vector<Element*> RTSolo::getLTRVec() const {
    std::vector<Element *> r = {ltr};
    return r;
}

const std::set<RT*> RTSolo::getNestSet() const {
    std::cerr << "Unsupported operation: A solo LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

const std::set<RT*> RTSolo::getOuterSet() const {
    return outerSet;
}

int RTSolo::getPPTStart() const {
    std::cerr << "Unsupported operation: A solo LTR RT does not have a PPT." << std::endl;
    throw std::exception();
}

int RTSolo::getPPTEnd() const {
    std::cerr << "Unsupported operation: A solo LTR RT does not have a PPT." << std::endl;
    throw std::exception();
}


std::pair<int, int> RTSolo::getLeftTSD() const {
    std::cerr << "Unsupported operation: A solo LTR RT does not have a TSD." << std::endl;
    throw std::exception();
}

std::pair<int, int> RTSolo::getRightTSD() const {
    std::cerr << "Unsupported operation: A solo LTR RT does not have a TSD." << std::endl;
    throw std::exception();
}

double RTSolo::getIdentityScore() const {
    std::cerr << "Unsupported operation: A solo LTR does not have an identity score." << std::endl;
    throw std::exception();
}

void RTSolo::setIsRC(bool isRC) {
    std::cerr << "Unsupported operation: Unable to set an LTR to be RC." << std::endl;
    throw std::exception();
}

void RTSolo::setPPT(int start, int end) {
    std::cerr << "Unsupported operation: A solo LTR RT does not have a PPT." << std::endl;
    throw std::exception();
}

void RTSolo::setTSD(int leftStart, int leftEnd, int rightStart, int rightEnd) {
    std::cerr << "Unsupported operation: A solo LTR RT does not have a TSD." << std::endl;
    throw std::exception();
}

void RTSolo::setIdentityScore(double identityScore) {
    std::cerr << "Unsupported operation: Unable to set a solo LTR to have an identity score." << std::endl;
    throw std::exception();
}

bool RTSolo::hasRightLTR() const
{
    return false;
}

bool RTSolo::hasLeftLTR() const
{
    return true;
}

void RTSolo::nest(RT* rt) {
    std::cerr << "Unsupported operation: A solo LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

void RTSolo::removeNest(RT* rt) {
    std::cerr << "Unsupported operation: A solo LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}



bool RTSolo::hasNest() const {
    std::cerr << "Unsupported operation: A solo LTR RT can not have nested elements." << std::endl;
    throw std::exception();
}

bool RTSolo::isNested() const {
    return outerSet.size() > 0 ? true : false;
}

bool RTSolo::couldNest(RT *rt) const {
    return false;
}

void RTSolo::addOuter(RT *rt){
    outerSet.insert(rt);
}

void RTSolo::removeOuter(RT *rt){
    if (outerSet.count(rt) == 1) {
        outerSet.erase(rt);
    }
}



void RTSolo::extend(int k, bool isForward) {
    if (isForward) {
        ltr->extendEnd(k);
    }
    else {
        ltr->extendStart(k);
    }
}

void RTSolo::expand(int pos, int length) {
    assert (pos >= 0);
    assert (length >= 0);
    if (pos <= ltr->getEnd()) {
        ltr->setEnd(ltr->getEnd() + length);
    }
    if (pos < ltr->getStart()) {
        ltr->setStart(ltr->getStart() + length);
    }

}

void RTSolo::push(int p) {

    int start = ltr->getStart();
    int end = ltr->getEnd();

    ltr->setEnd(end + p);
    ltr->setStart(start + p);

}

std::ostream &operator<<(std::ostream &os, const RTSolo &r)
{
    os << r.getStart() << " " << r.getEnd();
    return os;
}