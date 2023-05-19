/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Red
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

#include <string>
#include <filesystem>
#include <assert.h>
#include <optional>
#include <cmath>
#include <iostream>
#include <vector>

#include "../utility/Util.h"
#include "../utility/ILocation.h"
#include "../utility/Location.h"
#include "../nonltr/ChromListMaker.h"
#include "../nonltr/Trainer.h"
#include "../nonltr/Scanner.h"
#include "../nonltr/ChromosomeOneDigit.h"

class Red
{
private:
    // Variables
    std::string gnm;
    int cor;
    int k;
    int ord;
    double gau;
    double thr;
    int min;

    Trainer *trainer;


    // Methods
    void calcGnmLen();
    void calcMarkovOrd();
    void calcGauWidth();

public:
    // Constructor

    Red(std::string _gnm, std::optional<int> _cor = {}, std::optional<int> _k = {},
        std::optional<int> _ord = {}, std::optional<double> _gau = {},
        std::optional<double> _thr = {}, std::optional<int> _min = {},
        std::optional<std::string> _cnd = {});
    ~Red();

    // Getter and Setters
    int getK();

    // Methods
    /**
     * Print k-mer table to a file
     */
    void printTable(std::string fileName);

    /**
     * Print the trained HMM to a file
     */
    void printModel(std::string fileName);

    /**
     * Scan the genome to print the repeat locations, scores, and/or a fasta file with repeats masked
     * An optional dirPath may be given to print the chromosomes within also.
     * repeatpath, scorePath, and maskPath are output locations; dirPath is an input location
     */
    void scan(std::string repeatPath, int format, std::optional<std::string> scorePath = {}, 
                std::optional<std::string> maskPath = {}, std::optional<std::string> dirPath = {});


    /**
     * Score a seq (not too long)
     */
    std::vector<int> score(std::string &seq);

    /**
     * Percentage of count of non-zero scores over total count of scores
     * Tells the repetitiveness of a sequence  
     */
    double calcPercent(std::string &seq);

    /**
     * Calculates the median score of a sequence without zeros 
     */
    int calcMedianScore(std::string &seq);

    /**
     * Calculates the mean score of a sequence
     * If false is passed to withZero, then the mean is calculated without including 0's
    */
    double calcMeanScore(std::string &seq, bool withZero = true);


    const std::vector<ILocation*>* predictRepeats(std::string &seq);

};