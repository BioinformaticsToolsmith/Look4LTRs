/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * LocalAlignment
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

#include "LocalAlignment.h"

LocalAlignment::LocalAlignment(std::string &_seq1, std::string &_seq2, int _match, int _mismatch, int _gapOpen, int _gapContinue) : seq1(_seq1), seq2(_seq2)
{
    // seq1 is the first sequence
    // seq2 is the sequence to be aligned to seq1
    // match is the score for a match
    // mismatch is the score for a mismatch
    // gapOpen is the score for opening a gap
    // gapContinue is the score for continuing a gap

    match = _match;
    mismatch = _mismatch;
    gapOpen = _gapOpen;
    gapContinue = _gapContinue;

    // Initializing values for alignment
    score = 0;
    length = 0;
    sameCount = 0;
    alignEnd1 = 0;
    alignEnd2 = 0;

    assert(seq1.size() > 0);
    assert(seq2.size() > 0);
    assert(std::max({match, mismatch, gapOpen, gapContinue}) > 0);

    align();

}

LocalAlignment::~LocalAlignment()
{

}

void LocalAlignment::align() {


    int rows = seq1.size() + 1;
    int cols = seq2.size() + 1;

    // Building matrices on heap
    int** scoreMatrix = new int*[rows];
    int** delGapMatrix = new int*[rows];
    int** insGapMatrix = new int*[rows];

    for (int i = 0; i < rows; ++i) {
        scoreMatrix[i] = new int[cols]();
        delGapMatrix[i] = new int[cols]();
        insGapMatrix[i] = new int[cols]();
    }

    // Filling matrices
    for (int i = 1; i < seq1.size() + 1; i++) {
        for (int j = 1; j < seq2.size() + 1; j++) {

            // If the two bases are the same, then match, otherwise mismatch
            int m = seq1[i - 1] == seq2[j - 1] ? match : mismatch;
            m += scoreMatrix[i - 1][j - 1];


            // Maximum score of the match/mismatch, gaps, and 0
            delGapMatrix[i][j] = std::max({scoreMatrix[i - 1][j] + gapOpen, delGapMatrix[i - 1][j] + gapContinue});
            insGapMatrix[i][j] = std::max({scoreMatrix[i][j - 1] + gapOpen, insGapMatrix[i][j - 1] + gapContinue});
            scoreMatrix[i][j] = std::max({m, delGapMatrix[i][j], insGapMatrix[i][j], 0});

            // Updating the end position of the alignment if the score is higher than the max score
            if (scoreMatrix[i][j] > score) {
                score = scoreMatrix[i][j];
                alignEnd1 = i;
                alignEnd2 = j;
            }
        }
    }

    // Checking for out-of-bounds indexing
    assert(alignEnd1 < seq1.size() + 1);
    assert(alignEnd2 < seq2.size() + 1);
    
    alignStart1 = alignEnd1;
    alignStart2 = alignEnd2;
    
    while (scoreMatrix[alignStart1][alignStart2] > 0) {
        bool same = seq1[alignStart1 - 1] == seq2[alignStart2 - 1];
        if (scoreMatrix[alignStart1][alignStart2] == scoreMatrix[alignStart1 - 1][alignStart2 - 1] + (same ? match : mismatch) ) {
            // If the two bases are the same, then add one to sameCount for the similarity calculation
            if (same) {
                sameCount++;
            }
            alignStart1--;
            alignStart2--;
        }
        else if (scoreMatrix[alignStart1][alignStart2] == delGapMatrix[alignStart1][alignStart2]) {
            alignStart1--;
        }
        else {
            alignStart2--;
        }
        length++;
    }



    // freeing matrices on heap from memory
    for (int i = 0; i < rows; ++i) {
        delete[] scoreMatrix[i];
        delete[] delGapMatrix[i];
        delete[] insGapMatrix[i];
    }

    delete[] scoreMatrix;
    delete[] delGapMatrix;
    delete[] insGapMatrix;
    
}

int LocalAlignment::getScore() {
    return score;
}

int LocalAlignment::getLength() {
    return length;
}

int LocalAlignment::getSameCount() {
    return sameCount;
}

double LocalAlignment::getSimilarity() {
    return double(sameCount) / length;
}

std::pair<int, int> LocalAlignment::getAlignLoc1() {
    return std::make_pair(alignStart1, alignEnd1);
}

std::pair<int, int> LocalAlignment::getAlignLoc2() {
    return std::make_pair(alignStart2, alignEnd2);
}