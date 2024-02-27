/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Detector.h
 * 
 *  Created on: Oct 13, 2022
 *      Author: Anthony B. Garza.
 *    Reviewer: Hani Z. Girgis
 * 
 * Purpose: Detects if consecutive stretches should be merged, merges them into an element if so
 *          and returns a vector of these elements.  Uses StretchFeature to build a feature table
 *          and Stochastic Gradient Descent to predict
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

#include "SGD.h"
#include "Stretch.h"
#include "LtrUtility.h"
#include "Element.h"
#include "LtrParameters.h"
#include "../Matrix.h"
#include "LtrParameters.h"
#include "../red/Red.h"
#include "../utility/ILocation.h"

#include "assert.h"
#include <string>
#include <unordered_map>

class Detector {
    public:

        /**
         * Constructors 
         */
        Detector(Red &_red, std::string &_seq);
        Detector(Red &_red, std::string &_seq, std::string &otherSeq);
        ~Detector();

        /**
         * Methods
        */
        // Returns a vector of elements merged according to the SGD classifier
        std::vector<Element> apply(std::vector<Stretch> &stretchVec);

        // prints out the Element variables
        void printElements(std::vector<Element> &elementVec, int showAmount = -1);

        std::vector<int> getPrediction();

        Matrix getFeatures();

        Matrix getScaledFeatures();

        const std::vector<ILocation*>* getRepeats();

        Matrix buildMatrix(std::vector<Stretch> &stretchVec);

    private:

        /**
         * Variables
         */
        Red &red;

        std::string &seq;
        
        SGD classifier;
    
        const int featureCount = 10;
        
        // Mean and std vectors for standardization
        std::vector<double> &meanVec;
        std::vector<double> &stdVec;

        std::vector<int> scoreVec;
        const std::vector<utility::ILocation *> *locationVec;

        std::vector<int> prediction;
        Matrix features;
        Matrix fMatrix;

        // stretch -> the overlapping Red's segment
        // This stretch and the one after it overlap with Red's segment
        std::unordered_map<Stretch*, ILocation*> overlapMap;

        /**
         * Methods
        */
        Matrix scale(Matrix &m);

        std::vector<Element> mergeStretches(std::vector<Stretch> &stretchVec, std::vector<int> &predVec);

};