#pragma once

#include "RT.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <string>

#include <iostream>

class LtrUtility {
public:

    inline static bool isEqual(double x, double y){
        return fabs(x - y) < std::numeric_limits<double>::epsilon();
    }

    inline static bool isGreaterEqual(double x, double y) {
        return (std::fabs(x - y) < std::numeric_limits<double>::epsilon() || x >= y);
    }

    inline static bool isGreater(double x, double y) {
        return (std::fabs(x - y) > std::numeric_limits<double>::epsilon());
    }

    inline static bool isLess(double x, double y) {
        return (std::fabs(y - x) > std::numeric_limits<double>::epsilon());
    }


    // Calculates the median without zero's
    static double calcMedian(std::vector<int> &scoreVec);

    // Calculates the mean
    static double calcMean(std::vector<int> &scoreVec);

    // Calculates the percentage of non-zero scores over all scores
    static double calcPercent(std::vector<int> &scoreVec);

    // Given a string, returns the reverse complement
    static std::string reverseComplement(std::string &seq);

    /**
     * Given a string of DNA and two characters, any time the first character appears
     * in the string, replace it with the second. Return this new string.
     * This should be done as efficient memory-wise and time-wise as possible.
    */
    static std::string replaceChar(std::string &seq, char oldChar, char newChar);

    static int getAlignmentLength(std::string &seq1, std::string &seq2, int match, int mismatch, int gapOpen, int gapContinue);

    static void sortRTs(std::vector<RT*> &rtVec);

    static void rankRTs(std::vector<RT*> &rtVec, bool sort = true);

    static void removeDuplicateRTs(std::vector<RT*> &rtVec, bool sort = true, double overlap = 0.5);

    static void removeNests(std::vector<RT*> &rtVec);

    static void removeIllFormat(std::vector<RT*> &rtVec);

    static void renameCaseType(std::vector<RT*> &rtVec, std::string currName, std::string newName);

private:
    // A table of complements for each nucleotide
    static const std::unordered_map<char, char> complementTable;

};
