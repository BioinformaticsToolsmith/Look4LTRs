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
 *          and returns a vector of these elements.
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

#include "Detector.h"

/**
 * _red: is the Red API
 * _seq: is a chromosome sequence
 */
Detector::Detector(Red &_red, std::string &_seq)
    : red(_red), seq(_seq), meanVec(LtrParameters::MEAN_VECTOR), stdVec(LtrParameters::STD_VECTOR), classifier(LtrParameters::WEIGHT_VECTOR)
{
    scoreVec = red.score(seq);
    locationVec = red.predictRepeats(seq);
}

Detector::Detector(Red &_red, std::string &_seq, std::string &otherSeq)
    : red(_red), seq(_seq), meanVec(LtrParameters::MEAN_VECTOR), stdVec(LtrParameters::STD_VECTOR), classifier(LtrParameters::WEIGHT_VECTOR)
{
    scoreVec = red.score(otherSeq); 
    locationVec = red.predictRepeats(otherSeq); 
}

// [OK]
Detector::~Detector()
{
    for (int i = 0; i < locationVec->size(); i++)
    {
        delete locationVec->at(i);
    }
    delete locationVec;
}

/**
 * The classifier is applied twice: once on the forward stretches and once on the backward stretches.
 */
std::vector<Element> Detector::apply(std::vector<Stretch> &stretchVec)
{
    // This is the result
    std::vector<Element> elementVec;

    if (stretchVec.size() == 0)
    {
        // Do nothing
    }
    else if (stretchVec.size() == 1)
    {
        elementVec.push_back(Element(stretchVec[0]));
    }
    else
    {
        overlapMap.clear();
        features = buildMatrix(stretchVec);
        fMatrix = scale(features);
        prediction = classifier.predict(fMatrix);

        // Post condtion
        assert(prediction.size() == fMatrix.getNumRow());


        elementVec = mergeStretches(stretchVec, prediction);
    }
    int min = LtrParameters::MIN_LTR;
    elementVec.erase(std::remove_if(elementVec.begin(), elementVec.end(),
                                  [min](const Element& ele) {
                                      return ele.getSize() < min;
                                  }),
                   elementVec.end());

    elementVec.shrink_to_fit();

    return elementVec;
}

std::vector<Element> Detector::mergeStretches(std::vector<Stretch> &stretchVec, std::vector<int> &predVec) {
    std::vector<Element> r;
    bool previous = false;
    for (int i = 0; i < stretchVec.size() - 1; i++) {
        if (!previous) {
            r.push_back(Element{stretchVec.at(i)});
        }

        if (predVec[i] == 1) {
            r.back().merge(stretchVec.at(i + 1));
            previous = true;
        }
        else {
            previous = false;
        }
    }
    if (!previous) {
        r.push_back(Element{stretchVec.back()});
    }
    return r;
}


/**
 * Extract features from stretches
 */
Matrix Detector::buildMatrix(std::vector<Stretch> &stretchVec)
{
    // The matrix size is the number of stretches minus 1 (because features are
    // extracted from consecutive pairs) times the features.

    /**
     * There are 10 features:
     * First Stretch Size     : size of the first stretch
     * Second Stretch Size    : size of the second stretch
     * Gap Size               : size of the gap between the stretches
     * Height Difference      : difference of stretches' heights
     * Median Score Difference: difference of median (does not count zeros) Red scores of both stretches
     * Mean Score Difference  : difference of mean (involves zero in calculations) Red scores of both stretches
     * First Mean Score       : mean Red Score of the first stretch
     * Second Mean Score      : mean Red Score of the second stretch
     * Gap Mean Score         : mean Red Score of the gap between the stretches
     * Repeat overlap         : 1 if both stretches overlap the same repeat detected by Red
     */

    overlapMap.clear();

    // Feature matrix
    Matrix r;

    if (!stretchVec.empty())
    {
        r = Matrix(stretchVec.size() - 1, featureCount);

        // Extract the repeat-overlap feature
        int i = 0; // the strech index
        int j = 0; // the repeat index
        int label = -1;
        while (i < stretchVec.size() - 1)
        {
            if (j < locationVec->size())
            {
                if (stretchVec[i].getEnd() < locationVec->at(j)->getStart())
                {
                    label = 0;
                }
                else if (stretchVec[i + 1].getStart() > locationVec->at(j)->getEnd() + 1)
                {
                    j++;
                    continue;
                }
                else if (stretchVec[i].isPartOf(locationVec->at(j)->getStart(), locationVec->at(j)->getEnd() + 1, LtrParameters::MIN_STRETCH_OVERLAP) &&
                         stretchVec[i + 1].isPartOf(locationVec->at(j)->getStart(), locationVec->at(j)->getEnd() + 1, LtrParameters::MIN_STRETCH_OVERLAP))
                {
                    label = 1;
                }
                else
                {
                    label = 0;
                }
            }
            else
            {
                label = 0;
            }

            overlapMap[&stretchVec[i]] = label == 1 ? locationVec->at(j) : nullptr;
            i++;
        }

        // Extract the rest of the features
        for (int i = 0; i < stretchVec.size() - 1; i++)
        {
            r(i, 0) = stretchVec[i].getSize();
            r(i, 1) = stretchVec[i + 1].getSize();
            r(i, 2) = stretchVec[i].calculateGap(stretchVec[i + 1]);
            r(i, 3) = std::abs(stretchVec[i].getMedianHeight() - stretchVec[i + 1].getMedianHeight());

            // Extract Red scores
            std::vector<int> firstStretch(scoreVec.begin() + stretchVec[i].getStart(), scoreVec.begin() + stretchVec[i].getEnd());
            std::vector<int> secondStretch(scoreVec.begin() + stretchVec[i + 1].getStart(), scoreVec.begin() + stretchVec[i + 1].getEnd());
            std::vector<int> gapStretch(scoreVec.begin() + stretchVec[i].getEnd(), scoreVec.begin() + stretchVec[i + 1].getStart());

            double firstMean = LtrUtility::calcMean(firstStretch);
            double secondMean = LtrUtility::calcMean(secondStretch);
            double gapMean = gapStretch.size() > 0 ? LtrUtility::calcMean(gapStretch) : 0.0;

            r(i, 4) = std::abs(LtrUtility::calcMedian(firstStretch) - LtrUtility::calcMedian(secondStretch));
            r(i, 5) = std::abs(firstMean - secondMean);
            r(i, 6) = firstMean;
            r(i, 7) = secondMean;
            r(i, 8) = gapMean;
            r(i, 9) = overlapMap[&stretchVec.at(i)] != nullptr ? 1 : 0;
        }
    }

    return r;
}

/**
 * Standardize the data
 */
Matrix Detector::scale(Matrix &m)
{
    int rowCount = m.getNumRow();
    int colCount = m.getNumCol();

    Matrix r(rowCount, colCount); // result

    // All features but the Repeat Overlap
    for (int row = 0; row < rowCount; row++)
    {
        for (int col = 0; col < colCount - 1; col++)
        {
            r(row, col) = (m(row, col) - meanVec[col]) / stdVec[col];
        }
        r(row, featureCount - 1) = m(row, featureCount - 1);
    }

    return r;
}

/**
 * showAmount: the number of elements to print
 * if it is 0 or less, then print all elements. 
 * This method waits on the user input to print more elements.
 */
void Detector::printElements(std::vector<Element> &elementVec, int showAmount)
{
    if (showAmount <= 0)
    {
        for (auto ele : elementVec)
        {
            std::cout << ele << std::endl;
        }
    }
    else
    {
        int counter = 0;
        for (auto ele : elementVec)
        {
            std::cout << ele << std::endl;
            counter++;
            if (counter == showAmount)
            {
                counter = 0;
                std::cin.get();
            }
        }
    }
}

std::vector<int> Detector::getPrediction()
{
    return prediction;
}

Matrix Detector::getFeatures()
{
    return features;
}

Matrix Detector::getScaledFeatures() {
    return fMatrix;
}

const std::vector<ILocation *> *Detector::getRepeats()
{
    return locationVec;
}