/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * LocAlign
 * 
 *  Created on: X X, 20XX
 *      Author: Anthony B. Garza.
 * Reviewer:
 *   Purpose: Local Alignment class; allows user to query for information between two strings
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

#include <string>
#include <assert.h>
#include <algorithm>
#include <utility>
#include <iostream>

class LocalAlignment
{
private:
    // Variables

    // Constructor variables
    std::string &seq1;
    std::string &seq2;
    int match;
    int mismatch;
    int gapOpen;
    int gapContinue;
    bool findSequence;

    // Output variables
    int score;
    int length;
    int sameCount;
    int alignStart1;
    int alignStart2;
    int alignEnd1;
    int alignEnd2;

    // Methods
    void align();

public:
    
    // Constructor

    LocalAlignment(std::string &_seq1, std::string &_seq2, int match, int mismatch, int gapOpen, int gapContinue);
    ~LocalAlignment();

    // Getter and Setters
    int getScore();

    int getLength();

    int getSameCount();

    double getSimilarity();

    std::pair<int, int> getAlignLoc1();

    std::pair<int, int> getAlignLoc2();

    // Methods
};