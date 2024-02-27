#include "LtrUtility.h"

double LtrUtility::calcMedian(std::vector<int> &scoreVec) {
    double r;
    std::vector<int> noZeroVec;
    for (auto x : scoreVec) {
        if (x != 0) {
            noZeroVec.push_back(x);
        }
    }

    if (noZeroVec.size() > 0) {
        std::sort(noZeroVec.begin(), noZeroVec.end());
        int middle = noZeroVec.size() / 2;
        if (noZeroVec.size() % 2 == 1) {
            r = noZeroVec.at(middle);
        }
        else {
            r = (noZeroVec.at(middle) + noZeroVec.at(middle - 1)) / 2.0;
        }
    }
    else {
        r = 0.0;
    }

    return r;
}

double LtrUtility::calcMean(std::vector<int> &scoreVec) {
    double r;
    double sum = 0;
    for (auto x : scoreVec) {
        sum += x;
    }

    r = sum / scoreVec.size();
    return r;
}

double LtrUtility::calcPercent(std::vector<int> &scoreVec) {
    double count = 0.0;
    for (auto x : scoreVec) {
        if (x != 0) {
            count++;
        }
    }

    return count / scoreVec.size();

}

// Defining the complement table
const std::unordered_map<char, char> LtrUtility::complementTable = {
    {'A', 'T'},
    {'T', 'A'},
    {'C', 'G'},
    {'G', 'C'},
    {'N', 'N'}
};

// Given a string of DNA, returns the reverse complement
static std::string LtrUtility::reverseComplement(std::string &seq) {
    std::string r = seq;
    std::reverse(r.begin(), r.end());
    std::transform(r.begin(), r.end(), r.begin(), [](char c) {
        return complementTable.at(c);
    });
    return r;
}

/**
 * Given a string of DNA and two characters, any time the first character appears
 * in the string, replace it with the second. Return this new string.
*/
static std::string LtrUtility::replaceChar(std::string &seq, char oldChar, char newChar) {
    std::string r = seq;
    std::transform(seq.begin(), seq.end(), r.begin(), [oldChar, newChar](char c) {
        return c == oldChar ? newChar : c;
    });
    return r;
}

int LtrUtility::getAlignmentLength(std::string &seq1, std::string &seq2, int match, int mismatch, int gapOpen, int gapContinue) {
    int scoreMatrix[seq1.size() + 1][seq1.size() + 1] = {};
    int delGapMatrix[seq1.size() + 1][seq1.size() + 1] = {};
    int insGapMatrix[seq1.size() + 1][seq1.size() + 1] = {};

    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    for (int i = 1; i < seq1.size() + 1; i++) {
        for (int j = 1; j < seq2.size() + 1; j++) {
            int m = seq1[i - 1] == seq2[j - 1] ? match : mismatch;
            m += scoreMatrix[i - 1][j - 1];

            delGapMatrix[i][j] = std::max({scoreMatrix[i - 1][j] + gapOpen, delGapMatrix[i - 1][j] + gapContinue});
            insGapMatrix[i][j] = std::max({scoreMatrix[i][j - 1] + gapOpen, insGapMatrix[i][j - 1] + gapContinue});
            scoreMatrix[i][j] = std::max({m, delGapMatrix[i][j], insGapMatrix[i][j], 0});

            if (scoreMatrix[i][j] > max_score) {
                max_score = scoreMatrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }


    int i = max_i;
    int j = max_j;
    int locAlign = 0;

    while (scoreMatrix[i][j] > 0) {
        if (scoreMatrix[i][j] == scoreMatrix[i-1][j-1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch) ) {
            i--;
            j--;
        }
        else if (scoreMatrix[i][j] == delGapMatrix[i][j]) {
            i--;
        }
        else {
            j--;
        }
        locAlign++;
    }

    return locAlign;

}

static void LtrUtility::sortRTs(std::vector<RT*> &rtVec) {
    std::sort(rtVec.begin(), rtVec.end(), [](RT *r1, RT *r2)
    {
        bool order;
        if (r1->getStart() == r2->getStart()) {
            if (r1->getEnd() == r2->getEnd()) {
                order = r1->getCaseRank() > r2->getCaseRank();
            } else {
                order = r1->getEnd() < r2->getEnd();
            }
        } 
        else {
            order = r1->getStart() < r2->getStart();
        }
        return order;
    });
}


static void LtrUtility::rankRTs(std::vector<RT*> &rtVec, bool sort) {
    if (sort) {
        LtrUtility::sortRTs(rtVec);
    }
    for (int i = 0; i < rtVec.size(); i++) {
        auto &curr = rtVec.at(i);
        if (curr == nullptr) {
            continue;
        }

        for (int j = i + 1; j < rtVec.size(); j++) {
            // std::cout << "NEXT!" << std::endl;
            auto &next = rtVec.at(j);
            if (next == nullptr) {
                continue;
            }

            if (curr->getStart() == next->getStart() && curr->getEnd() == next->getEnd()) {
                if (curr->getCaseRank() >= next->getCaseRank()) {
                    delete next;
                    next = nullptr;
                }
            }
            else if (!curr->hasRightLTR() && next->hasRightLTR() && curr->calcOverlap(next) > 0) {
                delete curr;
                curr = nullptr;
                break;
            }
            else if (!next->hasRightLTR() && curr->hasRightLTR() && next->calcOverlap(curr) > 0) {
                delete next;
                next = nullptr;
            }
            else {
                break;
            }

        }
    }
    rtVec.erase(std::remove(rtVec.begin(), rtVec.end(), nullptr), rtVec.end());
}

static void LtrUtility::removeDuplicateRTs(std::vector<RT*> &rtVec, bool sort, double overlapThreshold) {
    if (sort) {
        LtrUtility::sortRTs(rtVec);
    }
    for (int i = 0; i < rtVec.size(); i++) {
        auto &curr = rtVec.at(i);
        if (curr == nullptr) {
            continue;
        }

        for (int j = i + 1; j < rtVec.size(); j++) {
            auto &next = rtVec.at(j);
            if (next == nullptr) {
                continue;
            }
            else if (next->getStart() >= curr->getEnd()) {
                break;
            }
            // If curr and next are RT Completes
            else if (curr->hasRightLTR() && next->hasRightLTR()) {
                // Get the overlap of the LTRs
                auto currLeft = curr->getLeftLTR();
                auto currRight = curr->getRightLTR();
                auto nextLeft = next->getLeftLTR();
                auto nextRight = next->getRightLTR();

                int leftOverlap = currLeft->calcOverlap(*nextLeft);
                int rightOverlap = currRight->calcOverlap(*nextRight);

                bool leftOverlapCond = isGreaterEqual(double(leftOverlap) / currLeft->getSize(), overlapThreshold) && isGreaterEqual(double(leftOverlap) / nextLeft->getSize(),overlapThreshold);
                bool rightOverlapCond = isGreaterEqual(double(rightOverlap) / currRight->getSize(), overlapThreshold) && isGreaterEqual(double(rightOverlap) / nextRight->getSize(),overlapThreshold);

                if (leftOverlapCond && rightOverlapCond) {
                    if (curr->getCaseRank() >= next->getCaseRank()) {
                        delete next;
                        next = nullptr;
                    }
                    else {
                        delete curr;
                        curr = nullptr;
                        break;
                    }
                }
            }
            // If curr and next are RT Solos
            else if (!curr->hasRightLTR() && !next->hasRightLTR()) {
                // Get the overlap of the LTRs
                int overlap = curr->calcOverlap(next);
                if (isGreaterEqual(double(overlap) / curr->getSize(), overlapThreshold) && isGreaterEqual(double(overlap) / next->getSize(),overlapThreshold)) {
                    if (curr->getCaseRank() >= next->getCaseRank()) {
                        delete next;
                        next = nullptr;
                    }
                    else {
                        delete curr;
                        curr = nullptr;
                        break;
                    }

                }
            }



        }
    }
    rtVec.erase(std::remove(rtVec.begin(), rtVec.end(), nullptr), rtVec.end());

}

static void LtrUtility::removeNests(std::vector<RT*> &rtVec) {
    for (auto &rt : rtVec) {
        if (rt->hasRightLTR() && rt->hasNest()) {
            for (auto &nest : rt->getNestSet()) {
                rt->removeNest(nest);
                nest->removeOuter(rt);
            }
        }
    }
}

static void LtrUtility::removeIllFormat(std::vector<RT*> &rtVec) {
    for (auto &rt : rtVec) {
        if (rt->hasRightLTR()) {
            int leftStart = rt->getLeftLTR()->getStart();
            int leftEnd = rt->getLeftLTR()->getEnd();
            int rightStart = rt->getRightLTR()->getStart();
            int rightEnd = rt->getRightLTR()->getEnd();
            if (leftStart > leftEnd || leftEnd > rightStart || rightStart > rightEnd) {
                delete rt;
                rt = nullptr;
            }
        }
    }
    rtVec.erase(std::remove(rtVec.begin(), rtVec.end(), nullptr), rtVec.end());
}


static void LtrUtility::renameCaseType(std::vector<RT*> &rtVec, std::string currName, std::string newName) {
    for (auto &rt : rtVec) {
        if (rt->getCaseType() == currName) {
            rt->setCaseType(newName);
        }
    }
}

static std::string LtrUtility::lcs(std::string &str1, std::string &str2) {
    int len1 = str1.size();
    int len2 = str2.size();
    
    // Create a matrix to store the lengths of subproblems
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1));
    
    int maxLength = 0; // To store length of the longest common substring
    int endIndex = 0; // To store the ending index of longest common substring in str1

    // Building the dp matrix
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            // If character at i-th index in str1 is equal to the character at j-th index in str2
            if (str1[i - 1] == str2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                
                // Update the maximum length and ending index
                if (dp[i][j] > maxLength) {
                    maxLength = dp[i][j];
                    endIndex = i - 1;
                }
            }
            // If character at i-th index in str1 is not equal to character at j-th index in str2
            else {
                dp[i][j] = 0;
            }
        }
    }
    
    // Return the longest common substring
    return str1.substr(endIndex - maxLength + 1, maxLength);
}


static std::string LtrUtility::getFileName(std::string filePath) {
    std::string name = filePath.substr(filePath.find_last_of('/') + 1);
    name = name.substr(0, name.find_last_of('.'));    
    return name;
}

static void LtrUtility::collectFastaFiles(std::vector<std::string> &collectVec, std::string &directory) {
    for (auto& p : std::filesystem::directory_iterator(directory)) {
        if (p.path().extension() == ".fasta" || p.path().extension() == ".fa") {
            collectVec.push_back(p.path().string());
        }
    }
}


static void LtrUtility::collectFastaFiles(std::vector<std::string> &collectVec, std::vector<std::string> &directoryVec) {
    for (auto& fastaDir : directoryVec) {
        collectFastaFiles(collectVec, fastaDir);
    }
}

static void LtrUtility::collectBedFiles(std::vector<std::string> &collectVec, std::string &directory) {
    for (auto& p : std::filesystem::directory_iterator(directory)) {
        if (p.path().extension() == ".bed") {
            collectVec.push_back(p.path().string());
        }
    }
}

static void LtrUtility::collectBedFiles(std::vector<std::string> &collectVec, std::vector<std::string> &directoryVec) {
    for (auto& bedDir : directoryVec) {
        collectBedFiles(collectVec, bedDir);
    }
}

static std::unique_ptr<IdentityCalculator<int32_t>> LtrUtility::buildCalculator(double threshold, std::string dbPath, int coreCount, bool skip) {
    SynDataGenerator dg{dbPath, threshold, coreCount};
    int64_t maxLength = dg.getMaxLength();
    assert (maxLength <= std::numeric_limits<int32_t>::max());

    std::unique_ptr<IdentityCalculator<int32_t>> ic = std::make_unique<IdentityCalculator<int32_t>>(&dg, coreCount, threshold, skip, true);
    return ic;
}
