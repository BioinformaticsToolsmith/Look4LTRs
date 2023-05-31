/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * PostProcess
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

#include "PostProcess.h"

PostProcess::PostProcess(std::vector<RT *> &_rtVec, std::vector<RT *> &_complexVec, int _max, IdentityCalculator<int32_t> &_ic, std::string *_seq) : rtVec(_rtVec), complexVec(_complexVec), ic(_ic)
{
    max = _max;
    seq = _seq;
}

PostProcess::~PostProcess()
{

}

void PostProcess::apply(bool nestOnly) {
    nest();
    if (!nestOnly) {
        extendK(rtVec);
        extendK(complexVec);
        extendIdentity(rtVec);
    }
}

void PostProcess::nest() {
    // Top of the stack is the last, most outer RT to attempt to nest within.
    std::stack<RT*> nestStack;

    for (int i = 0; i < rtVec.size(); i++) {
        // If nestStack is empty, there is no current outer RT to check
        while (!nestStack.empty()) {
            // Nest the RT at rtVec[i] in the outer RT at nestStack's top IF it is within the interior.
            if (nestStack.top()->couldNest(rtVec[i])) {
                nestStack.top()->nest(rtVec[i]);
                break;
            }
            // If the RT at rtVec[i] is after the end of the RT at nestStack's top, pop the top of the stack. No more RTs can be logically nested within it.
            else if (rtVec[i]->getStart() >= nestStack.top()->getEnd() || rtVec[i]->getEnd() >= nestStack.top()->getEnd()) {
                nestStack.pop();
            }
            else {
                break;
            }
        }
        // Now we push the current RT onto the stack to check if the next RTs can be nested within it.
        nestStack.push(rtVec[i]);
    }
    // for (int i = 0; i < rtVec.size(); i++) {
    //     for (int j = i + 1; j < rtVec.size(); j++) {
    //         if (rtVec[i]->couldNest(rtVec[j])) {
    //             rtVec[i]->nest(rtVec[j]);
    //         }
    //         else if (rtVec[j]->getStart() >= rtVec[i]->getEnd()) {
    //             break;
    //         }
    //     }
    // }
}

void PostProcess::extendLTR(RT *rtPtr) {
    
    Element *leftLTR = rtPtr->getLeftLTR();
    Element *rightLTR = rtPtr->getRightLTR();

    // Get the smaller of the two elements and the larger of the two elements;
    // store it in a smallPtr and largePtr
    Element *smallPtr = leftLTR;
    Element *largePtr = rightLTR;
    if (leftLTR->getSize() > rightLTR->getSize()) {
        smallPtr = rightLTR;
        largePtr = leftLTR;
    }

    // Get the new starts of the smaller and larger LTRs
    auto positions = largePtr->getHypoStarts(smallPtr);
    if (positions.first != std::numeric_limits<int>::max()) {
        assert (positions.first >= largePtr->getStart() && positions.first < largePtr->getEnd());
        assert (positions.second >= smallPtr->getStart() && positions.second < smallPtr->getEnd());

        // Based on the new start of the larger, how much is left in the LTR?
        int largeRemainderForward = largePtr->getEnd() - positions.first;
        int largeRemainderBackward = positions.first - largePtr->getStart() + 1;

        // Based on the new start of the smaller, where should it now end?
        int smallEnd = std::max(positions.second + largeRemainderForward, smallPtr->getEnd());
        if (smallEnd > seq->size()) {
            smallEnd = seq->size();
        }
        int smallStart = std::min(positions.second - largeRemainderBackward, smallPtr->getStart());
        if (smallStart < 0) {
            smallStart = 0;
        }

        double originalScore = 0.0;
        std::string largeLTRSeq = seq->substr(largePtr->getStart(), largePtr->getSize());

        // if (!ic.isImpossible(largePtr->getSize(), smallPtr->getSize(), LtrParameters::MIN_IDENTITY)) {
        //     std::string smallLTRSeq = seq->substr(smallPtr->getStart(), smallPtr->getSize());
        //     originalScore = ic.score(&largeLTRSeq, &smallLTRSeq);
        // }
        // // Get new small LTR's sequence
        // std::string smallLTRSeq = seq->substr(smallStart, smallEnd - smallStart);

        std::string smallLTRSeq = seq->substr(smallPtr->getStart(), smallPtr->getSize());
        originalScore = smallLTRSeq.size() >= ic.getK()? ic.score(&largeLTRSeq, &smallLTRSeq) : 0.0;
        smallLTRSeq = seq->substr(smallStart, smallEnd - smallStart);
        
        double extendedScore = ic.score(&largeLTRSeq, &smallLTRSeq);

        // Compare with Identity
        if (extendedScore > originalScore) {
        // if (!LtrUtility::isEqual(extendedScore, 0.0) && extendedScore > originalScore) {
            smallPtr->setStart(smallStart);
            smallPtr->setEnd(smallEnd);
            rtPtr->setIdentityScore(extendedScore);
        }
        else {
            rtPtr->setIdentityScore(originalScore);
        }
    }
    else {
        std::string largeLTRSeq = seq->substr(largePtr->getStart(), largePtr->getSize());
        std::string smallLTRSeq = seq->substr(smallPtr->getStart(), smallPtr->getSize());
        double identityScore = smallLTRSeq.size() >= ic.getK()? ic.score(&largeLTRSeq, &smallLTRSeq) : 0.0;
        rtPtr->setIdentityScore(identityScore);
    }

}

void PostProcess::extendIdentity(std::vector<RT*> &vec) {

    for (auto rtPtr : vec) {
        if (rtPtr->hasRightLTR()) {
            extendLTR(rtPtr);
        }
    }





    // for (auto rtPtr : vec) {
    //     if (rtPtr->hasRightLTR()) {
    //         Element *left = rtPtr->getLeftLTR();
    //         Element *right = rtPtr->getRightLTR();

    //         // Regions in the rightLTR that have a match in the leftLTR
    //         auto rightRegion = left->findMissingMatch(right);

    //         // Extending the right LTR backwards
    //         if (rightRegion.first != -1) {

    //             // Asking if the new start doesn't overlap with the left LTR
    //             if (rightRegion.first > left->getEnd()) {
    //                 std::string leftLTRSeq = seq->substr(left->getStart(), left->getSize());
    //                 std::string newRightLTRSeq = seq->substr(rightRegion.first, right->getEnd() - rightRegion.first);

    //                 if (!LtrUtility::isEqual(ic.score(&leftLTRSeq, &newRightLTRSeq), 0.0)) {
    //                     right->extendStart(right->getStart() - rightRegion.first);
    //                 }
    //             }
    //         }

    //         // Extending the right LTR forwards
    //         if (rightRegion.second != -1) {
                
    //             // Asking if the new end doesn't extend past the end of the sequence;
    //             if (rightRegion.second <= seq->size()) {
    //                 std::string leftLTRSeq = seq->substr(left->getStart(), left->getSize());
    //                 std::string newRightLTRSeq = seq->substr(right->getStart(), rightRegion.second - right->getStart());

    //                 if (!LtrUtility::isEqual(ic.score(&leftLTRSeq, &newRightLTRSeq), 0.0)) {
    //                     right->extendEnd(rightRegion.second - right->getEnd());
    //                 }        
    //             }
    //         }

    //         // continue;
    //         // Extending the left LTR backwards
    //         auto leftRegion = right->findMissingMatch(left);
    //         if (leftRegion.first != -1) {
                
    //             // Asking if the new start is above 0
    //             if (leftRegion.first >= 0) {

    //                 std::string rightLTRSeq = seq->substr(right->getStart(), right->getSize());
    //                 std::string newLeftLTRSeq = seq->substr(leftRegion.first, left->getEnd() - leftRegion.first);
    //                 if (!LtrUtility::isEqual(ic.score(&rightLTRSeq, &newLeftLTRSeq), 0.0)) {
    //                     left->extendStart(left->getStart() - leftRegion.first);
    //                 }
    //             }
    //         }
    //         // Extending the left LTR forwards
    //         if (leftRegion.second != -1) {
                
    //             // Asking if the new end doesn't overlap with the right LTR
    //             if (leftRegion.second < right->getStart()) {

    //                 std::string rightLTRSeq = seq->substr(right->getStart(), right->getSize());
    //                 std::string newLeftLTRSeq = seq->substr(left->getStart(), leftRegion.second - left->getStart());
    //                 if (!LtrUtility::isEqual(ic.score(&rightLTRSeq, &newLeftLTRSeq), 0.0)) {
    //                     left->extendEnd(leftRegion.second - left->getEnd());
    //                 }        
    //             }
    //         }
    //     }
    // }
}

void PostProcess::extendK(std::vector<RT*> &vec) {
    // @@@@@ ACCOUNT FOR SEQ SIZE
    for (auto rtPtr : vec) {
        if (rtPtr->getEnd() + LtrParameters::K - 1 <= max) {
            rtPtr->extend(LtrParameters::K - 1);
        } 
        else {
            rtPtr->extend(max - rtPtr->getEnd());
        }
        // rtPtr->extend(LtrParameters::K - 1);

    }
}