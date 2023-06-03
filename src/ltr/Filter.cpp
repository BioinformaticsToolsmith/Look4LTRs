/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Filter
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

#pragma once

#include "Filter.h"

Filter::Filter(std::vector<RT*> &_rtVec, Red &_red, std::string *_seq) : rtVec(_rtVec), red(_red), seq(_seq)
{
}

void Filter::apply() {
    // Check to see if rtVec has any RTs to filter
    if (rtVec.size() > 0) {
        filterLength();
        filterN();
        filterLengthRatio();
        // filterExtendLength();
        filterPPT();
        filterMITE();
        filterLonesomes();

        markTSD();

    }
}

void Filter::filterLength() {
    // vector to hold the RTs that pass the filter
    std::vector<RT* > r;
    for (auto rtPtr : rtVec) {

        bool isFail = true;

        int leftLtrLength = rtPtr->getLeftLTR()->getSize();
        bool isLeftLtrFit = leftLtrLength >= LtrParameters::MIN_LTR && leftLtrLength <= LtrParameters::MAX_LTR;
        bool isRightLtrFit = true;
        bool isInteriorFit = true;

        if (rtPtr->hasRightLTR()) {

            int interiorLength = rtPtr->getRightLTR()->getStart() - rtPtr->getLeftLTR()->getEnd();
            int trueInteriorLength;

            auto nestSet = rtPtr->getNestSet();


            if (nestSet.empty()) {
                trueInteriorLength = interiorLength;
            }
            else {
                std::vector<RT*> nestVec = std::vector<RT*>(nestSet.begin(), nestSet.end());
                std::sort(nestVec.begin(), nestVec.end(), [](RT* a, RT* b) { return a->getLeftLTR()->getStart() < b->getLeftLTR()->getStart(); });

                std::vector<std::pair<int, int>> nestRanges;
                for (auto nestPtr : nestVec) {
                    if (nestPtr->hasRightLTR()) {
                        nestRanges.push_back(std::make_pair(nestPtr->getLeftLTR()->getStart(), nestPtr->getRightLTR()->getEnd()));
                    }
                    else {
                        nestRanges.push_back(std::make_pair(nestPtr->getLeftLTR()->getStart(), nestPtr->getLeftLTR()->getEnd()));
                    }
                }

                std::vector<std::pair<int, int>> mergedRanges;
                for (auto &j : nestRanges) {
                    if (mergedRanges.empty()) {
                        mergedRanges.push_back(j);
                    }
                    else {
                        if (j.first <= mergedRanges.back().second) {
                            mergedRanges.back().second = std::max(mergedRanges.back().second, j.second);
                        }
                        else {
                            mergedRanges.push_back(j);
                        }
                    }
                }

                int mergedLength = 0;
                for (auto &j : mergedRanges) {
                    mergedLength += j.second - j.first;
                }

                trueInteriorLength = interiorLength - mergedLength;
            }

            int rightLtrLength = rtPtr->getRightLTR()->getSize();

            rtPtr->getNestSet();

            isInteriorFit = interiorLength >= LtrParameters::MIN_INTERIOR;
            isRightLtrFit = rightLtrLength >= LtrParameters::MIN_LTR && rightLtrLength <= LtrParameters::MAX_LTR;
        }

        if (isLeftLtrFit && isRightLtrFit && isInteriorFit) {
            isFail = false;
        }

        // if (isLeftLtrFit) {
        //     if (rtPtr->hasRightLTR()) {
        //         /**
        //          * The length of the RT should exclude the length of any nested elements for the filter.
        //          * getRange() will get the ranges of an LTR RT excluding any nested element.
        //          */
        //         int rtLength = 0;
                // for (auto &j : rtPtr->getRange()) {
        //             rtLength += j.second - j.first;
        //         }
        //         int rightLtrLength = rtPtr->hasRightLTR()? rtPtr->getRightLTR()->getSize() : 1.0;
        //         int interiorLength = rtLength - leftLtrLength - rightLtrLength;

        //         int trueInteriorLength = rtPtr->getRightLTR()->getStart() - rtPtr->getLeftLTR()->getEnd();

        //         bool isRtFit = rtLength >= LtrParameters::MIN_RT && rtLength <= LtrParameters::MAX_RT;
        //         bool isRightLtrFit = rightLtrLength >= LtrParameters::MIN_LTR && rightLtrLength <= LtrParameters::MAX_LTR;
        //         bool isInteriorFit = interiorLength >= LtrParameters::MIN_INTERIOR && interiorLength <= LtrParameters::MAX_INTERIOR && trueInteriorLength >= LtrParameters::MIN_INTERIOR && trueInteriorLength <= LtrParameters::MAX_INTERIOR;

        //         if (isRtFit && isRightLtrFit && isInteriorFit) {
        //             isFail = false;
        //         }
        //     }
        //     else {
        //         if (rtPtr->getSize() >= LtrParameters::MIN_LTR) {
        //             isFail = false;
        //         }
        //     }

        // }



        if (isFail) {
            removeNests(rtPtr);
            delete rtPtr;
        }
        else {
            r.push_back(rtPtr);
        }
    }

    rtVec = r;
}

void Filter::filterN() {
    std::vector<RT* > r; // vector to hold the RTs that pass the filter
    // int nCount = 0;
    for (auto &rt : rtVec) {
        bool leftLtrPass = false;
        bool rightLtrPass = false;

        std::string leftLtrSeq = seq->substr(rt->getLeftLTR()->getStart(), rt->getLeftLTR()->getSize());
        if (std::count(leftLtrSeq.begin(), leftLtrSeq.end(), 'N') / double(rt->getLeftLTR()->getSize()) <= LtrParameters::MAX_N_RATIO) {
            leftLtrPass = true;
        }

        if (rt->hasRightLTR()) {
            std::string rightLtrSeq = seq->substr(rt->getRightLTR()->getStart(), rt->getRightLTR()->getSize());
            if (std::count(rightLtrSeq.begin(), rightLtrSeq.end(), 'N') / double(rt->getRightLTR()->getSize()) <= LtrParameters::MAX_N_RATIO) {
                rightLtrPass = true;
            }
        }
        else {
            rightLtrPass = true;
        }

        if (leftLtrPass && rightLtrPass) {
            r.push_back(rt);
        }
        else {
            // nCount++;
            removeNests(rt);
            delete rt;
            rt = nullptr;
        }
    }
    // std::cout << "Removed " << nCount << " RTs due to N ratio." << std::endl;
    rtVec = r;
}

void Filter::filterIdentity() {
    std::vector<RT* > r; // vector to hold the RTs that pass the filter
    for (auto &rt : rtVec) {
        if (rt->hasRightLTR()) {
            if (rt->getIdentityScore() >= LtrParameters::MIN_IDENTITY_FILTER) {
                r.push_back(rt);
            }
            else {
                removeNests(rt);
                delete rt;
                rt = nullptr;
            }
        }
        else {
            r.push_back(rt);
        }
    }
    rtVec = r;
}

void Filter::filterLengthRatio() {
    std::vector<RT* > r; // vector to hold the RTs that pass the filter
    for (auto &rt : rtVec) {
        if (rt->hasRightLTR()) {
            double l1 = rt->getLeftLTR()->getSize();
            double l2 = rt->getRightLTR()->getSize();
            double ratio = l1 < l2 ? l1 / l2 : l2 / l1;
            if (ratio >= LtrParameters::MIN_LENGTH_RATIO) {
                r.push_back(rt);
            }
            else {
                int leftStart = rt->getLeftLTR()->getStart();
                int leftEnd = rt->getLeftLTR()->getEnd();
                int rightStart = rt->getRightLTR()->getStart();
                int rightEnd = rt->getRightLTR()->getEnd();

                std::string leftLtrSeq = seq->substr(leftStart, leftEnd - leftStart);
                std::string rightLtrSeq = seq->substr(rightStart, rightEnd - rightStart);
                LocalAlignment la(leftLtrSeq, rightLtrSeq, 2, -3, -5, -2);
                int length = la.getLength();
                double similarity = la.getSimilarity();

                if (length >= LtrParameters::MIN_LTR && similarity >= LtrParameters::MIN_IDENTITY) {
                    int newLeftStart, newLeftEnd, newRightStart, newRightEnd;
                    std::tie(newLeftStart, newLeftEnd) = la.getAlignLoc1();
                    std::tie(newRightStart, newRightEnd) = la.getAlignLoc2();


                    rt->getLeftLTR()->setStart(leftStart + newLeftStart);
                    rt->getLeftLTR()->setEnd(leftStart + newLeftEnd);
                    rt->getRightLTR()->setStart(rightStart + newRightStart);
                    rt->getRightLTR()->setEnd(rightStart + newRightEnd);
                    r.push_back(rt);

                }
                else {
                    removeNests(rt);
                    delete rt;
                    rt = nullptr;
                }
            }
        }
        else {
            r.push_back(rt);
        }
    }

    rtVec = r;
}

// void Filter::filterLengthRatio() {
//     std::vector<RT* > r; // vector to hold the RTs that pass the filter
//     for (auto &rt : rtVec) {
//         if (rt->hasRightLTR()) {
//             int leftStart = rt->getLeftLTR()->getStart();
//             int leftEnd = rt->getLeftLTR()->getEnd();
//             int rightStart = rt->getRightLTR()->getStart();
//             int rightEnd = rt->getRightLTR()->getEnd();

//             std::string leftLtrSeq = seq->substr(leftStart, leftEnd - leftStart);
//             std::string rightLtrSeq = seq->substr(rightStart, rightEnd - rightStart);
//             LocalAlignment la(leftLtrSeq, rightLtrSeq, 2, -3, -5, -2);
//             int length = la.getLength();
//             double similarity = la.getSimilarity();

//             int newLeftStart, newLeftEnd, newRightStart, newRightEnd;


//             bool pass = false;
//             if (length >= LtrParameters::MIN_LTR && similarity >= LtrParameters::MIN_IDENTITY) {
//                 std::tie(newLeftStart, newLeftEnd) = la.getAlignLoc1();
//                 std::tie(newRightStart, newRightEnd) = la.getAlignLoc2();
//                 newLeftStart = newLeftStart + leftStart;
//                 newLeftEnd = newLeftEnd + leftStart;
//                 newRightStart = newRightStart + rightStart;
//                 newRightEnd = newRightEnd + rightStart;

//                 leftLtrSeq = seq->substr(newLeftStart, newLeftEnd - newLeftStart);
//                 rightLtrSeq = seq->substr(newRightStart, newRightEnd - newRightStart);

//                 // Count N's
//                 if (std::count(leftLtrSeq.begin(), leftLtrSeq.end(), 'N') / double(leftLtrSeq.size()) <= LtrParameters::MAX_N_RATIO &&
//                     std::count(rightLtrSeq.begin(), rightLtrSeq.end(), 'N') / double(rightLtrSeq.size()) <= LtrParameters::MAX_N_RATIO) {
//                     pass = true;
//                 }

//             }
//             if (pass) {
//                 rt->getLeftLTR()->setStart(newLeftStart);
//                 rt->getLeftLTR()->setEnd(newLeftEnd);
//                 rt->getRightLTR()->setStart(newRightStart);
//                 rt->getRightLTR()->setEnd(newRightEnd);
//                 r.push_back(rt);
//             }
//             else {
//                 removeNests(rt);
//                 delete rt;
//                 rt = nullptr;
//             }
//         }
//         else {
//             r.push_back(rt);
//         }
//     }

//     rtVec = r;
// }

// void Filter::filterExtendLength() {
//     std::vector<RT* > r; // vector to hold the RTs that pass the filter
//     int lengthCount = 0;
//     int alignmentCount = 0;
//     int failCount = 0;
//     for (auto &rt : rtVec) {
//         if (rt->hasRightLTR()) {
//             double l1 = rt->getLeftLTR()->getSize();
//             double l2 = rt->getRightLTR()->getSize();
//             double ratio = l1 < l2 ? l1 / l2 : l2 / l1;
//             if (ratio >= LtrParameters::MIN_LENGTH_RATIO) {
//                 lengthCount++;
//                 r.push_back(rt);
//             }
//             else {
//                 int leftStart = rt->getLeftLTR()->getStart();
//                 int leftEnd = rt->getLeftLTR()->getEnd();
//                 int rightStart = rt->getRightLTR()->getStart();
//                 int rightEnd = rt->getRightLTR()->getEnd();
//                 bool foundMatch = false;

//                 int sizeDifference = std::abs(l1 - l2);
//                 Element *left = rt->getLeftLTR();
//                 Element *right = rt->getRightLTR();
//                 Element *smaller = l1 < l2 ? left : right;
//                 Element *larger = l1 < l2 ? right : left;

//                 double leftSimilarity = 0.0;
//                 int leftInteriorLength = -1;
//                 double rightSimilarity = 0.0;
//                 int rightInteriorLength = -1;

//                 // Extend smaller left
//                 int newSmallStart = smaller->getStart() - sizeDifference;
//                 if (!(newSmallStart <= leftEnd && smaller == right)) {
//                     std::string smallLtrSeq = seq->substr(newSmallStart, smaller->getSize() + sizeDifference);
//                     std::string largeLtrSeq = seq->substr(larger->getStart(), larger->getSize());
//                     LocalAlignment la(smallLtrSeq, largeLtrSeq, 2, -3, -5, -2);
//                     leftInteriorLength = smaller == right ? newSmallStart - leftEnd : -1;
//                     if (leftInteriorLength >= LtrParameters::MIN_INTERIOR) {
//                         leftSimilarity = la.getSimilarity();
//                         foundMatch = true;
//                     }
//                 }

//                 // Extend smaller right
//                 int newSmallEnd = smaller->getEnd() + sizeDifference;
//                 if (!(newSmallEnd >= rightStart && smaller == left)) {
//                     std::string smallLtrSeq = seq->substr(smaller->getStart(), smaller->getSize() + sizeDifference);
//                     std::string largeLtrSeq = seq->substr(larger->getStart(), larger->getSize());
//                     LocalAlignment la(smallLtrSeq, largeLtrSeq, 2, -3, -5, -2);
//                     rightInteriorLength = smaller == left ? rightStart - newSmallEnd : -1;
//                     if (rightInteriorLength >= LtrParameters::MIN_INTERIOR) {
//                         rightSimilarity = la.getSimilarity();
//                         foundMatch = true;
//                     }
//                 }

//                 if (foundMatch && (leftSimilarity >= LtrParameters::MIN_IDENTITY || rightSimilarity >= LtrParameters::MIN_IDENTITY)) {
//                     if (leftSimilarity > rightSimilarity) {
//                         smaller->setStart(newSmallStart);
//                     }
//                     else {
//                         smaller->setEnd(newSmallEnd);
//                     }
//                     alignmentCount++;
//                     r.push_back(rt);
//                 }
//                 else {
//                     failCount++;
//                     removeNests(rt);
//                     delete rt;
//                     rt = nullptr;
//                 }
//             }
//         }
//         else {
//             r.push_back(rt);
//         }
//     }

//     rtVec = r;
    

//     std::cout << "Length: " << lengthCount << " Alignment: " << alignmentCount << " Fail: " << failCount << std::endl;
// }

// void Filter::filterLengthRatio() {
//     std::vector<RT* > r; // vector to hold the RTs that pass the filter
//     for (auto &rt : rtVec) {
//         if (rt->hasRightLTR()) {
//             double l1 = rt->getLeftLTR()->getSize();
//             double l2 = rt->getRightLTR()->getSize();
//             double ratio = l1 < l2 ? l1 / l2 : l2 / l1;
//             if (ratio >= LtrParameters::MIN_LENGTH_RATIO) {
//                 r.push_back(rt);
//             }
//             else {
//                 removeNests(rt);
//                 delete rt;
//                 rt = nullptr;
//             }
//         }
//         else {
//             r.push_back(rt);
//         }
//     }
//     rtVec = r;
// }


void Filter::filterRepetivity() {
    std::vector<RT* > r;
    for (auto rtPtr : rtVec) {
        // Checking if the whole RT without nested elements is repetitive
        // Only if the RT is not a solo LTR check:
        // Left LTR repetivity
        // Right LTR repetivity
        // Interior (no nested elements) repetivitiy
        bool isFail = true;
        auto rtScoreVec = scoreSeq(rtPtr->getStart(), rtPtr->getEnd());

        auto s1 = rtScoreVec.begin();
        auto s2 = s1 + rtPtr->getLeftLTR()->getSize();
        

        std::vector<int> leftLtrScoreVec{s1, s2};
        bool isLeftLtrRep = LtrUtility::calcPercent(leftLtrScoreVec) >= LtrParameters::MIN_PERC;

        if (isLeftLtrRep) {
            if (rtPtr->hasRightLTR()) {

                auto s3 = s2 + rtPtr->getRightLTR()->getStart() - rtPtr->getLeftLTR()->getEnd();
                
                std::vector<int> rightLtrScoreVec{s3, rtScoreVec.end()};
                bool isRightLtrRep = LtrUtility::calcPercent(rightLtrScoreVec) >= LtrParameters::MIN_PERC;

                // First and last index contains the LTR ranges
                Range rangeVec = rtPtr->getRange();
                assert(rangeVec.size() >= 3);

                int start = rtPtr->getStart();
                std::vector<int> interiorScoreVec;
                for (int i = 1; i < rangeVec.size() - 1; i++) {
                    auto locs = rangeVec.at(i);
                    interiorScoreVec.insert(interiorScoreVec.end(), rtScoreVec.begin() + locs.first - start , rtScoreVec.begin() + locs.second - start);
                }
                bool isInteriorRep =interiorScoreVec.size() > LtrParameters::MIN_INTERIOR ? LtrUtility::calcPercent(interiorScoreVec) >= LtrParameters::MIN_PERC : false;

                if (isRightLtrRep && isInteriorRep) {
                    isFail = false;
                }

            }
            else {
                isFail = false;
            }
        }

        if (isFail) {
            removeNests(rtPtr);
            delete rtPtr;
        }
        else {
            r.push_back(rtPtr);
        }
    }

    rtVec = r;
}

/**
 * This method will determine if a given sequence is a MITE.
 * The sequence is a MITE if the following conditions are met:
 * 1. The sequence is between 2 * MIN_MITE_TIR_SIZE and MAX_MITE_SIZE (2 Minimum TIR size because a MITE must be able to hold two TIRs)
 * 2. The first TIR_SEARCH_RANGE nucleotides of the sequence align with the reverse complement of the last TIR_SEARCH_RANGE nucleotides.
 * 2a. This alignment must be at least MIN_MITE_TIR_SIZE nucleotides long.
 * 2b. Based on the alignment, there must be at least MITE_SIMILARITY similarity between the two TIRS.
 */
bool Filter::findMITE(std::string &ele) {
    bool r = false;

    // Minimum and maximum size; A MITE must be able to hold at least two TIRs
    if (ele.size() >= LtrParameters::MIN_MITE_TIR_SIZE * 2 && ele.size() <= LtrParameters::MAX_MITE_SIZE) {

        std::string leftTirSeq;
        std::string rightTirSeq;

        // Get up to TIR_SEARCH_RANGE from both sides if the element is big enough
        if (ele.size() > LtrParameters::TIR_SEARCH_RANGE * 2) {
            leftTirSeq = ele.substr(0, LtrParameters::TIR_SEARCH_RANGE);
            rightTirSeq = ele.substr(ele.size() - LtrParameters::TIR_SEARCH_RANGE, LtrParameters::TIR_SEARCH_RANGE);
        }
        // Otherwise, get up to half the element for both sides.
        else {
            int half = ele.size()/2;
            leftTirSeq = ele.substr(0, half);
            rightTirSeq = ele.substr(ele.size() - half, half);
        }

        // Reverse complement the right side; the right side of a TIR is the reverse complement of the left side
        std::string rightTirSeqRev = LtrUtility::reverseComplement(rightTirSeq);

        // Perform local alignment
        LocalAlignment la(leftTirSeq, rightTirSeqRev, 2, -3, -5, -2);
        int length = la.getLength();
        double similarity = la.getSimilarity();
        
        // If the alignment length is greater than LtrParameters::MIN_MITE_TIR_ALIGN and the similarity is greater than LtrParameters::MITE_SIMILARITY, then it is a MITE
        if (length >= LtrParameters::MIN_MITE_TIR_SIZE && LtrUtility::isGreaterEqual(similarity, LtrParameters::MITE_SIMILARITY)) {
            r = true;
        }
    }
    return r;
}

void Filter::filterMITE() {
    std::vector<RT* > r; // vector to hold the RTs that pass the filter

    for (auto rtPtr : rtVec) {

        // Only check for complete RTs
        if (rtPtr->hasRightLTR()) {

            // Getting the sequences of the left and right LTRs
            auto leftLTR = rtPtr->getLeftLTR();
            auto rightLTR = rtPtr->getRightLTR();
            std::string leftSeq = seq->substr(leftLTR->getStart(), leftLTR->getSize());
            std::string rightSeq = seq->substr(rightLTR->getStart(), rightLTR->getSize());

            // Checking for MITES in both; only drop if both are MITES (helps prevent dropping of true positives)
            if (findMITE(leftSeq) && findMITE(rightSeq)) {
                // remove any nests, both outer and inner
                removeNests(rtPtr);
                delete rtPtr;
                continue;
            }
        }
        r.push_back(rtPtr);
    }
    
    rtVec = r;
}

std::pair<int, int> Filter::findPPT(int border, int size, char replaceChar, char withChar) {
    std::pair<int, int> r = std::make_pair(-1, -1);

    bool isLarge = size > LtrParameters::MAX_PPT_DISTANCE * 2; // If we need to search half or up to max ppt distance
    int half = size / 2;

    int pptBorder = isLarge? std::abs(border + LtrParameters::MAX_PPT_DISTANCE) : std::abs(border + half);

    int smaller = std::min({std::abs(border), pptBorder});
    int larger = std::max({std::abs(border), pptBorder});
    std::string pptCandidate = seq->substr(smaller, larger - smaller);

    // Replacing every G from the candidate PPT region with A (C with T if reverse complement)
    std::string replaceSeq = LtrUtility::replaceChar(pptCandidate, replaceChar, withChar);
    // Creating a string as big as a PPT can be, but with all A's (C with T if reverse complement)
    std::string seqAll = std::string(LtrParameters::MAX_PPT_SIZE, withChar);

    // Getting local alignment
    LocalAlignment la(replaceSeq, seqAll, 2, -3, -5, -2);
    int length = la.getLength();

    // Only if the alignment length is at least the size of the PPT
    if (length >= LtrParameters::MIN_PPT_SIZE) {
        int pptStart, pptEnd;
        std::tie(pptStart, pptEnd) = la.getAlignLoc1();
        std::string ppt = pptCandidate.substr(pptStart, pptEnd - pptStart);

        int replaceCount = std::count(replaceSeq.begin(), replaceSeq.end(), replaceChar);
        int withCount = std::count(ppt.begin(), ppt.end(), withChar);

        if (withCount > replaceCount) {
            if (pptStart)
            r.first = smaller + pptStart;
            r.second = smaller + pptEnd;
        }
    }
    return r;
    
}


/**
 * Given each RT, search the interior for a poly purine tract (PPT)
 * This polypurine tract should be a bunch of purines (A or G) in a row
 * To search for this, we will take the interior, convert every G to an A,
 * create a sequence the size of the interior comprised of only A's,
 * Then perform alignment between the converted interior and the sequence of A's.
 * If the alignment is long enough, then we have a PPT.
 * The PPT can also be on the other side of the RT if the RT was reverse complemented.
 
 */
void Filter::filterPPT() {
    // vector of RTs that pass the filter
    std::vector<RT* > r;

    for (auto rtPtr : rtVec) {
        if (rtPtr->hasRightLTR()) {

            auto size = rtPtr->getIntSize();
            if (size < LtrParameters::MIN_PPT_SIZE) {
                removeNests(rtPtr);
                delete rtPtr;
                continue;
            }
            auto intStart = rtPtr->getLeftLTR()->getEnd();  // start of the interior
            auto intEnd = rtPtr->getRightLTR()->getStart(); // end of the interior

            auto rightPPT = findPPT(intEnd * -1, size, 'G', 'A');
            if (rightPPT.first != -1) {
                rtPtr->setPPT(rightPPT.first, rightPPT.second);
                r.push_back(rtPtr);
                continue;
            }

            auto leftPPT = findPPT(intStart, size, 'C', 'T');
            if (leftPPT.first != -1 && leftPPT.second <= intEnd) {
                rtPtr->setIsRC();
                rtPtr->setPPT(leftPPT.first, leftPPT.second);
                r.push_back(rtPtr);
                continue;
            }
            else {
                removeNests(rtPtr);
                delete rtPtr;
                continue;
            }

        }
        else {
            r.push_back(rtPtr);
        }
    }
    
    rtVec = r;

} 

/**
 * Find target site duplications. Due to complications with the merging of elements earlier,
 * we will check for TSDs within a range around the start of the 5' LTR and the end of the 3' LTR.
 * We will start by obtaining the sequence around the 5' LTR's start and the sequence around the 3' LTR's end.
 * Then we will perform local alignment.
 * Note that this method only marks rts with TSDs, it does not remove them.
*/
void Filter::markTSD() {

    for (auto rtPtr : rtVec) {

        // Only check if the RT is complete, i.e., has both LTRs
        if (rtPtr->hasRightLTR()) {
            // Get left TSD candidate's boundaries
            int leftTSDStart = std::max({0, rtPtr->getLeftLTR()->getStart() - LtrParameters::MAX_TSD_DISTANCE});
            int leftTSDEnd = rtPtr->getLeftLTR()->getStart();

            // Get right TSD candidate's boundaries
            int rightTSDStart = rtPtr->getRightLTR()->getEnd();
            int rightTSDEnd = std::min({int(seq->size()), rtPtr->getRightLTR()->getEnd() + LtrParameters::MAX_TSD_DISTANCE});

            // Getting the candidate TSD sequences
            std::string rightTSDSeq = seq->substr(rightTSDStart, rightTSDEnd - rightTSDStart);
            std::string leftTSDSeq = seq->substr(leftTSDStart, leftTSDEnd - leftTSDStart);

            if (rightTSDSeq.size() > 0 && leftTSDSeq.size() > 0) {

                LocalAlignment la(leftTSDSeq, rightTSDSeq, 2, -3, -5, -2);
                int length = la.getLength();

                // Only if the alignment length is at least the size of the smallest TSD.
                if (length >= LtrParameters::MIN_TSD_SIZE) {
                    int leftStart, leftEnd, rightStart, rightEnd;
                    std::tie(leftStart, leftEnd) = la.getAlignLoc1();
                    std::tie(rightStart, rightEnd) = la.getAlignLoc2();
                    leftStart += leftTSDStart;
                    leftEnd += leftTSDStart;
                    rightStart += rightTSDStart;
                    rightEnd += rightTSDStart;

                    rtPtr->setTSD(leftStart, leftEnd, rightStart, rightEnd);

                }
            }
        }
    }

}

void Filter::filterLonesomes() {
    std::vector<RT* > r;

    std::unordered_set<int> graphSet;

    for (auto rtPtr : rtVec) {
        if (rtPtr->hasRightLTR()) {
            graphSet.insert(rtPtr->getGraphGroup());
        }
    }
    for (auto rtPtr : rtVec) {
        if (rtPtr->hasRightLTR()) {
            r.push_back(rtPtr);
        }
        else {
            if (graphSet.count(rtPtr->getGraphGroup()) == 0) {
                removeNests(rtPtr);
                delete rtPtr;
            }
            else {
                r.push_back(rtPtr);
            }
        }
    }

    rtVec = r;
}

void Filter::removeNests(RT* rtPtr) {
    RT* outer = nullptr;
    if (rtPtr->hasOuter()) {
        outer = rtPtr->getOuterNest();
        outer->removeNest(rtPtr);
    }
    if (rtPtr->hasRightLTR()) {
        for (auto p : rtPtr->getNestSet())  {
            rtPtr->removeNest(p);
            if (outer != nullptr) {
                outer->nest(p);
            }
        }
    }
}

std::vector<int> Filter::scoreSeq(int s, int e) {
    std::string region = seq->substr(s, e - s);
    return red.score(region);
}


std::vector<RT*> Filter::getRtVec() {
    return rtVec;
}
