
#pragma once

#include "CaseSingle.h"

CaseSingle::CaseSingle(IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, Red &_red, const std::string *_seq) : CaseMatcher(_ic, _icRecent, _red, _seq)
{
    name = "Single";
    rank = 100;
}

void CaseSingle::apply(DirectedGraph<Element> &graph, std::vector<Element *> &forwardVec, std::vector<Element *> &backwardVec, int graphIndex)
{ /**
   * Simple Single RT
   *
   * >-->
   *
   * O
   *  \
   *   \
   *    O
   *
   * Simplest form of an RT.  A forward node matches a backward node AND the area between is repetitive throughout the genome.
   */

    try{
    for (auto fElePtr : forwardVec) {

        for (auto bElePtr : getDiagonals(fElePtr, graph)) {

            Element *newEle = nullptr;
            Element *leftLtr = nullptr;
            Element *rightLtr = nullptr;

            std::tie(newEle, leftLtr, rightLtr) = buildLtrs(fElePtr, bElePtr, graph);

            if (leftLtr != nullptr && rightLtr != nullptr && checkSingle(leftLtr, rightLtr)) {
                leftLtr = leftLtr == newEle? leftLtr : new Element{*leftLtr};
                rightLtr = rightLtr == newEle? rightLtr : new Element{*rightLtr};
                // std::tie(leftLtr, rightLtr) = assignLtrs(newEle, leftLtr, rightLtr);
                rtVec.push_back(new RTComplete{leftLtr, rightLtr, name, rank, graphIndex});
            }
            else if (newEle != nullptr) {
                delete newEle;
            }

        }
    }

    //filterSequential();
    } catch (...) {
        std::cout << "Single" << std::endl;
        throw std::exception();
    }
}

bool CaseSingle::checkSingle(Element* leftLtr, Element* rightLtr) {
    // Sequences for the LTRs and the interior

    bool r = false;
    bool isGap = leftLtr->calcGap(*rightLtr) > 0;

    if (isGap) {
        std::vector<int> rtScoreVec = scoreSeq(leftLtr->getStart(), rightLtr->getEnd());

        auto s1 = 0;
        auto s2 = s1 + leftLtr->getSize();
        auto s3 = s2 + rightLtr->getStart() - leftLtr->getEnd();


        std::vector<int> leftLtrScoreVec{rtScoreVec.begin() + s1, rtScoreVec.begin() + s2};
        std::vector<int> interiorScoreVec{rtScoreVec.begin() + s2, rtScoreVec.begin() + s3};
        std::vector<int> rightLtrScoreVec{rtScoreVec.begin() + s3, rtScoreVec.end()};
        r = checkSingle(leftLtrScoreVec, interiorScoreVec, rightLtrScoreVec);
    }


    // Any RT candidate has to have an interior, thus the gap between the two LTRs must not be 0


    return r;
}

bool CaseSingle::checkSingle(std::vector<int> &leftLtrScoreVec, std::vector<int> &interiorScoreVec, std::vector<int> &rightLtrScoreVec) {

    bool r = false;

    // Median calculation does not involve the zero RED scores
    double leftLtrScoreMedian = LtrUtility::calcMedian(leftLtrScoreVec);
    double interiorScoreMedian = LtrUtility::calcMedian(interiorScoreVec);
    double rightLtrScoreMedian = LtrUtility::calcMedian(rightLtrScoreVec); 

    double interiorPerc = LtrUtility::calcPercent(interiorScoreVec);

    r = interiorPerc >= LtrParameters::MIN_PERC;

    return r;
}

void CaseSingle::filterSequential() {

    if (rtVec.size() > 1) {

        std::vector<RT*> r;

        // RT -> score, where score is how many times the interior of the RT key matches the interior of another RT
        std::unordered_map<RT*, int> confidenceMap;
        for (auto rt : rtVec) {
            confidenceMap[rt] = 0;
        }

        for (int i = 0; i < rtVec.size() - 1; i++) {

            int iStart = rtVec.at(i)->getLeftLTR()->getEnd();
            int iEnd = rtVec.at(i)->getRightLTR()->getStart();
            std::string interior = seq->substr(iStart, iEnd - iStart);

            for (int j = i + 1; j < rtVec.size(); j++) {

                int cStart = rtVec.at(j)->getLeftLTR()->getEnd();
                int cEnd = rtVec.at(j)->getRightLTR()->getStart();
                std::string check = seq->substr(cStart, cEnd - cStart);

                bool isSame = calcSizeRatio(iStart, iEnd, cStart, cEnd) >= LtrParameters::MIN_IDENTITY? !LtrUtility::isEqual(ic.score(&interior, &check), 0.0) : false; 

                if (isSame){
                    confidenceMap[rtVec.at(i)]++;
                    confidenceMap[rtVec.at(j)]++;
                }
            }
        }

        // check for any confidence score above zero
        bool hasConfidence = false;
        for (auto rt : rtVec) {
            if (confidenceMap[rt] > 0) {
                hasConfidence = true;
                break;
            }
        }

        if (hasConfidence) {
            for (int i = 0; i < rtVec.size(); i++) {
                //std::cout << confidenceMap[rtVec.at(i)] << " ";
                if (confidenceMap[rtVec.at(i)] == 0) {
                    delete rtVec.at(i);
                    rtVec.at(i) = nullptr;
                }
                else {
                    r.push_back(rtVec.at(i));

                }
            }
            // std::cout << std::endl;
            rtVec = r;
        }

    }
}
