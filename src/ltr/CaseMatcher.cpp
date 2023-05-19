#include "CaseMatcher.h"

// [OK]
CaseMatcher::CaseMatcher(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq) : ic(_ic), icRecent(_icRecent), red(_red), seq(_seq)
{
}

// [OK]
bool CaseMatcher::isSame(Element *elementOne, Element *elementTwo, DirectedGraph<Element> &graph)
{
    bool r = false;

    if (graph.isBidirectional(*elementOne, *elementTwo) &&
        LtrUtility::isEqual(graph.retrieveWeight(*elementOne, *elementTwo), 0.0) &&
        LtrUtility::isEqual(graph.retrieveWeight(*elementTwo, *elementOne), 0.0))
    {
        r = true;
    }

    return r;
}

// []
// bool CaseMatcher::hasOverlap(Element *e, DirectedGraph<Element> &graph)
// {

//     bool r = false;
//     auto connectedVec = graph.retrieveConnectedValues(*e);
//     bool isVertical = false;
//     for (auto ptr : connectedVec)
//     {
//         if (e->isOverlap(*ptr) && e->getEnd() != ptr->getStart())
//         {
//             if (graph.retrieveWeight(*e, *ptr) != 0)
//             {
//                 r = true;
//                 break;
//             }
//             else if (!isVertical)
//             {
//                 isVertical = true;
//             }
//             else if (isVertical)
//             {
//                 r = true;
//                 break;
//             }
//         }
//     }

//     return r;
// }

// bool CaseMatcher::isOverlap(std::vector<Element *> &eVec, DirectedGraph<Element> &graph)
// {
//     bool r = false;
//     for (auto ptr : eVec)
//     {
//         if (hasOverlap(ptr, graph))
//         {
//             r = true;
//             break;
//         }
//     }

//     return r;
// }

// [OK]
// Enforces the first 80 of the 80-80-80 rule (length)
bool CaseMatcher::checkLength(std::vector<Element *> &eVec) {
    bool r = true;
    for (auto e : eVec) {
        if (!checkLength(e)) {
            r = false;
            break;
        }
    }
    return r;
}

// [OK]
bool CaseMatcher::checkLength(Element *e) {
    return e->getSize() >= LtrParameters::MIN_LTR;
}

// [OK]
int CaseMatcher::adjustSize(int s, int e) {
    int end = e + red.getK() - 1;
    if(end > seq->size()){
        end = seq->size();
    }
    return end - s;
}

// [OK]
std::vector<int> CaseMatcher::scoreSeq(int s, int e) {
    int size = adjustSize(s, e);
    std::string subSeq = seq->substr(s, size);
    return red.score(subSeq);
}

// [OK]
bool CaseMatcher::isAfter(Element *e1, Element *e2) {
    return e1->getStart() >= e2->getEnd();
}

// [OK]
Element * CaseMatcher::getVertical(Element *e, DirectedGraph<Element> &graph) {
    auto eleVec = graph.retrieveConnectedValues(*e);
    Element *r = nullptr;
    for (auto ele : eleVec) {
        if (LtrUtility::isEqual(graph.retrieveWeight(*e, *ele), 0.0)) {
            r = ele;
            break;
        } 
    }
    return r;
}

// std::vector<Element *> CaseMatcher::getTraversable(Element *e, DirectedGraph<Element> &graph) {
//     std::vector<Element *> r;
//     std::unordered_set<Element *> visitedSet;
//     std::stack<Element*> eleStack;
//     auto conVec = getDiagonals();
// }

// [OK]
/**
 * Merge the two vertical connections
 * @@@@@ We need to think about this step because it could be the source of the hyperextension problem
 */
Element * CaseMatcher::buildElement(Element *e, DirectedGraph<Element> &graph, bool &isDelete) {
    Element *vert = getVertical(e, graph);
    Element *r = e;
    if (vert != nullptr) {
        r = new Element{*e, *vert};
        isDelete = true;
    }
    return r;
}

std::tuple<Element*, Element*> CaseMatcher::assignLtrs(Element *built, Element *left, Element *right) {
    Element *l = nullptr;
    Element *r = nullptr;

    if (left == built) {
        l = built;
        r = new Element{*right};
    }
    else if(right == built) {
        l = new Element{*left};
        r = built;
    }
    else {
        throw std::exception();
    }

    return std::make_tuple(l, r);
}


/**
 * Check for recipricol coverage, if fails -> check for hyper merge and check recipricol coverage again.
 * Return three elements: 
 */
std::tuple<Element*, Element*, Element*> CaseMatcher::buildLtrs(Element *f, Element *b, DirectedGraph<Element> &graph) {
    // Finding the first RT's LTRs
    Element *newEle   = nullptr;
    Element *leftLtr  = nullptr;
    Element *rightLtr = nullptr;

    if (graph.isBidirectional(*f, *b)) {
        if (isReciprocal(f, b, graph)) {
            leftLtr = f;
            rightLtr = b;
        }
        // Only if there is a possibility of hyperextension.
        else {
            bool delLeft = false;
            bool delRight = false;
            newEle = isStretchReciprocal(f, b, graph, delLeft, delRight);
            if (newEle != nullptr) {
                assert(delLeft != delRight);

                leftLtr  = delLeft?  newEle : f;
                rightLtr = delRight? newEle : b;

                assert (leftLtr != rightLtr);
            }
        }
    }

    return std::make_tuple(newEle, leftLtr, rightLtr);
}


// [OK]
Element * CaseMatcher::retrieveFirst(Element *e1, Element *e2) {
    return e1->getStart() < e2->getStart() ? e1 : e2;
}

// For identity
// [OK]
double CaseMatcher::calcSizeRatio(int s1, int e1, int s2, int e2) {
    double size1 = e1 - s1;
    double size2 = e2 - s2;
    return size1 > size2? size2/size1 : size1/size2; 
}

// [OK]
double CaseMatcher::calcSizeRatio(int size1, int size2) {
    return size1 > size2? size2/double(size1) : size1/double(size2);
}

bool CaseMatcher::isReciprocal(Element *left, Element *right, DirectedGraph<Element>& graph) {
    assert (graph.isBidirectional(*left, *right));
    assert (left->getIsForward());
    assert (!right->getIsForward());

    double fWeight = graph.retrieveWeight(*left, *right);
    double bWeight = graph.retrieveWeight(*right, *left);

    return LtrUtility::isGreaterEqual(fWeight, LtrParameters::MIN_WEIGHT) && LtrUtility::isGreaterEqual(bWeight, LtrParameters::MIN_WEIGHT);
}

Element* CaseMatcher::isStretchReciprocal(Element *left, Element *right, DirectedGraph<Element>& graph, bool &delLeft, bool &delRight) {
    assert (graph.isBidirectional(*left, *right));
    assert (left->getIsForward());
    assert (!right->getIsForward());

    Element* r = nullptr;

    double fWeight = graph.retrieveWeight(*left, *right);
    double bWeight = graph.retrieveWeight(*right, *left);

    if (LtrUtility::isGreaterEqual(fWeight, LtrParameters::MIN_WEIGHT) || LtrUtility::isGreaterEqual(bWeight, LtrParameters::MIN_WEIGHT)) {

        Element *hyper = nullptr;
        Element *ele = nullptr;

        if (fWeight > bWeight) {
            hyper = right;
            ele = left;
            delRight = true;
        }
        else if (fWeight < bWeight){
            hyper = left;
            ele = right;
            delLeft = true;
        }
        
        if (ele != nullptr) {

            std::vector<Stretch*> stretchVec = hyper->getStretchVec();
            int start = -1;
            int end = -1;

            std::vector<Stretch*> keepStretchVec;

            for (auto stretch : stretchVec) {
                Stretch hypotheticalMatch = stretch->buildMatch();
                if (ele->calcOverlap(hypotheticalMatch) > 0) {
                    start = start == -1? hypotheticalMatch.getStart() : std::min(start, hypotheticalMatch.getStart());
                    end = end == -1? hypotheticalMatch.getEnd() : std::max(end, hypotheticalMatch.getEnd());
                    keepStretchVec.push_back(stretch);
                }
            }

            if (start >= 0) {
                assert (end > start);

                double coverage = ele->calcOverlap(start, end) / (end - start);

                if (LtrUtility::isGreaterEqual(coverage, LtrParameters::MIN_WEIGHT)) {
                    r = new Element(keepStretchVec);
                }
            }

        }
    }

    return r;
}

std::vector<Element*> CaseMatcher::getDiagonals(Element *e, DirectedGraph<Element> & graph) {
    auto conVec = graph.retrieveConnectedValues(*e);
    std::vector<Element *> r;
    for (auto ele : conVec) {
        if (!LtrUtility::isEqual(graph.retrieveWeight(*e, *ele), 0.0)) {
            r.push_back(ele);
        }
    }

    std::sort(r.begin(), r.end(), [](Element *ele1, Element *ele2) {return ele1->getStart() < ele2->getStart();});
    return r;

}

std::vector<Element*> CaseMatcher::getHorizontals(Element *e, DirectedGraph<Element> &graph, bool forward) {
    std::vector<Element*> r;

    // Go forward by position
    if ((e->getIsForward() == true && forward == true) || (e->getIsForward() == false && forward == false)) {
        auto conVec = getDiagonals(e, graph);
        for (auto ele : conVec) {
            auto h = getVertical(ele, graph);
            if (h != nullptr) {
                r.push_back(h);
            }
        }
    }
    // Go backward by position
    else {
        Element *v = getVertical(e, graph);
        if (v != nullptr) {
            for (auto ele : getDiagonals(v, graph)) {
                r.push_back(ele);
            }
        }

    }

    std::sort(r.begin(), r.end(), [](Element *ele1, Element *ele2) {return ele1->getStart() < ele2->getStart();});

    return r;
}

std::vector<Element*> CaseMatcher::getHorizontalsAcross(Element *e, DirectedGraph<Element> &graph, bool forward) {
    std::unordered_set<Element*> set;
    std::vector<Element*> r;
    auto horizontalVec = getHorizontals(e, graph, forward);
    set.insert(horizontalVec.begin(), horizontalVec.end());

    for (auto horiz : horizontalVec) {
        auto sideVec = getHorizontals(horiz, graph, forward);
        set.insert(sideVec.begin(), sideVec.end());
    }

    r.insert(r.end(), set.begin(), set.end());

    return r;
}


std::string CaseMatcher::substrInterior(Element * e1, Element *e2, std::string *seq) {
    assert (e1->getEnd() <= e2->getStart());
    return seq->substr(e1->getEnd(), e2->getStart() - e1->getEnd());
}

bool CaseMatcher::areSeqSame(std::string &seq1, std::string &seq2) {
    return calcSizeRatio(seq1.size(), seq2.size()) >= LtrParameters::MIN_IDENTITY && 
            !LtrUtility::isEqual(ic.score(&seq1, &seq2), 0.0);
}

bool CaseMatcher::areRecentSeqSame(std::string &seq1, std::string &seq2) {
    return calcSizeRatio(seq1.size(), seq2.size()) >= LtrParameters::MIN_IDENTITY && 
            !LtrUtility::isEqual(icRecent.score(&seq1, &seq2), 0.0);
}

std::vector<RT *> CaseMatcher::getRTVec() const
{
    return rtVec;
}

// bool CaseMatcher::hasRTs() const
// {
//     return rtVec.size() > 0 ? true : false;
// }

int CaseMatcher::getRank() const {
    return rank;
}

std::string CaseMatcher::getName() const {
    return name;
}
