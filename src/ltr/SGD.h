/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * SGD.h
 * 
 *  Created on: Oct 5, 2022
 *      Author: Anthony B. Garza
 *      Reviewer: Hani Z. Girgis, PhD
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

#include "../Matrix.h"

#include <vector>
#include <cstdlib>
#include <random>
#include <assert.h>

class SGD {
    public:

		// For training and setting biases and weights to random values
		// static const int MIN_RAND = 0;
		// static const int MAX_RAND = 100;


        /**
         * Constructor
         */
        SGD(std::vector<double> weightVec); // the first weight is the bias
        // Constructor for training
        // SGD(double learningRate, int iterationCount);

        /**
         * Methods
         */ 
        
        // Given a matrix of features and a vector of labels, fit
        // void fit(Matrix &featureMatrix, std::vector<bool> labelVec);

        // Given a matrix of features, predict
        std::vector<int> predict(Matrix &featureMatrix);

        /**
         * Getter and Setters
         */

        std::vector<double> getWeightVec();
        void setWeightVec(std::vector<double>& weightVec);

		//        int getFeatureCount();

    private:
        /**
         * Variables
         */ 


		//        // learning rate when training
	//        double learningRate;
	//
	//        // number of iterations when training
	//        int iterationCount = 100;
	//
	//        // y-intercept, starting bias for predictions
	//        double bias;

        /*
         * Weights of each feature
         * First index contains the bias
         */ 
        std::vector<double> weightVec;
        
        // Number of features; should be the same as the size of weights - 1
        int featureCount;
};
