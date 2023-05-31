#include "SGD.h"
#include "StretchFeatures.h"
#include "Stretch.h"
#include "../Matrix.h"
#include "LtrParameters.h"

#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <iostream>

Stretch makeStretch(int start, int end, int height, int mark = 0);
Matrix readFeatureFile(std::string featurePath);

int main() {

    std::string featurePath = "/home/transposons/Projects/Identity/cases/sgd/case3.csv";
    //std::string featurePath;
    //std::cin >> featurePath;


    StretchFeatures a(LtrParameters::MEAN_VECTOR, LtrParameters::STD_VECTOR);
    Matrix features = readFeatureFile(featurePath);
    Matrix scaledFeatures = a.scale(features);

    // Predicting with user-given  weights
    SGD classifier(LtrParameters::WEIGHT_VECTOR);

    std::vector<int> prediction = classifier.predict(scaledFeatures);

    for (auto x: prediction) {
        std::cout << x << std::endl;
    }

    

    return 0;
}

Stretch makeStretch(int start, int end, int height, int mark) {
    Stretch a = {start, end, mark};
    a.setMedianHeight(height);
    return a;
}

Matrix readFeatureFile(std::string featurePath) {
    std::ifstream featureFile(featurePath);
    
    std::vector<std::array<int,4>> featureVec;

    int colCount = 4;
    int rowCount = 0; // To be defined when reading;
    if (!featureFile.is_open()) {
        throw std::exception();
        std::cerr << featurePath << " could not be opened!";
    }
    std::string line;
    std::getline(featureFile, line);  // Skipping header

    while (std::getline(featureFile, line)) {
        rowCount++;  // Getting number of rows;
        std::istringstream stream(line);
        std::string val;
        std::array<int, 4> featureArray;
        std::getline(stream, val, ',');
        std::getline(stream, val, ',');


        std::getline(stream, val, ',');
        featureArray[0] = std::stoi(val);

        std::getline(stream, val, ',');
        featureArray[1] = std::stoi(val);

        std::getline(stream, val, ',');
        featureArray[2] = std::stoi(val);

        std::getline(stream, val);
        featureArray[3] = std::stoi(val);

        featureVec.push_back(featureArray);
    }

    Matrix r(rowCount, colCount);
    for (int row = 0; row < rowCount; row++) {
        for (int col = 0; col < colCount; col++) {
            r(row, col) = featureVec[row][col];
        }
    }
    
    return r;
}
