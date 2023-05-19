//#include "../utility/LocAlign.h"
#include <string>
#include <iostream>
#include <algorithm>
#include "../ltr/LocalAlignment.h"

int getLength(std::string &seq1, std::string &seq2, int match, int mismatch, int gapOpen, int gapContinue);

int main() {
    // std::string seq1 = "TGACGTCCGCAGTCCACCAGTACCATGCTGTATCCGCCCCAAGGTGATCTATCAAGTGTAGTGCTAATTGGGATTCGAGCGCCACGTACAACTCAAGAGA";
    // std::string seq2 = "ACATCTATTCAAATAGCGTCCGTACCTAGATCGCTTCGCAATTCGTCGGCCGCGCTGCCGTTAAGCGCTCAACTTATGCTCTAGTAGGGCCCTGACAACG";
    // std::string seq1 = "ATCG";
    // std::string seq2 = "ATCG";
    // std::string seq1 = "CGTAGGCTTAAGGTTA";
    // std::string seq2 = "ATAGATA";
    // std::string seq1 = "CACGTGATCAA";
    // std::string seq2 = "AGCATCGGTTG";

    // std::string seq = "gggtcttgttaacttgtgccctaagggcacaagataagaaactaaatgtagaaacttaatcttgaaaattgtgcattcataaccttaaaagcttaaaaacttatattttcaatgcaatatttattttttaattctattttaagatccttatcttgtgccctaagggcacaagttaacatttctc";
    // std::string seq1 = "gggtcttgttaacttgtgccctaagggcacaagataagaaactaaatgtagaaacttaatcttgaaaattgtgcattcataaccttaaaagc";
    // std::string seq2 = "gagaaatgttaacttgtgcccttagggcacaagataaggatcttaaaatagaattaaaaaataaatattgcattgaaaatataagtttttaa";
    // std::string seq = "   caggggtcggcaaactatggcccgagggccaaatctggcccgccgcctgtttttgtacggcccgcgagctaagaatggtttttacatttttaaatacaataaaactttatttnaaaatgtaaaaaccattcttagctcatgggccgtacaaaaacaggcggtaggccggatttggcctgcgcaggccgtagtttgccgacccctg";
    // std::string seq1 = "aactttatttnaaaatgtaaaaaccattcttagctcatgggccgtacaaaaacaggcggtaggccggatttggcctgcgcaggccgtagtttgccgacccctg";
    // std::string seq2 = "ttattgtatttaaaaatgtaaaaaccattcttagctcgcgggccgtacaaaaacaggcggcgggccagatttggccctcgggccatagtttgccgacccctg";
    // std::string seq1 = "CCCCAATAAGGG";
    // std::string seq2 = "AAAAAAAAAA";
    std::string seq1 = "AAGTAGGAAG";
    std::string seq2 = "AAATAAAAAA";

    LocalAlignment la(seq1, seq2, 1, -1, -2, -2);
    std::cout << la.getLength() << std::endl;
    std::cout << la.getScore() << std::endl;
    std::cout << la.getAlignLoc1().first << " " << la.getAlignLoc1().second << std::endl;
    std::cout << la.getAlignLoc2().first << " " << la.getAlignLoc2().second << std::endl;

    // LocAlign locAlign(seq1.c_str(), 0, seq1.size() - 1, seq2.c_str(), 0, seq2.size() - 1, 2, -3, 5, 2);
    // std::cout << locAlign.getLength() << std::endl;
    // locAlign.printAlignment();

    // std::cout << std::endl;

    // getLength(seq1, seq2, 2, -3, -5, -2);

}

int getLength(std::string &seq1, std::string &seq2, int match, int mismatch, int gapOpen, int gapContinue) {
    int scoreMatrix[seq1.size() + 1][seq2.size() + 1] = {};
    int delGapMatrix[seq1.size() + 1][seq2.size() + 1] = {};
    int insGapMatrix[seq1.size() + 1][seq2.size() + 1] = {};

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

    std::string align1;
    std::string align2;
    int i = max_i;
    int j = max_j;

    while (scoreMatrix[i][j] > 0) {
        if (scoreMatrix[i][j] == scoreMatrix[i-1][j-1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch) ) {
            align1 += seq1[i - 1];
            align2 += seq2[j - 1];
            i--;
            j--;
        }
        else if (scoreMatrix[i][j] == delGapMatrix[i][j]) {
            align1 += seq1[i - 1];
            align2 += "-";
            i--;
        }
        else {
            align1 += "-";
            align2 += seq2[j - 1];
            j--;
        }
    }

    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());

    std::cout << max_i << " " << i << std::endl;
    std::cout << max_j << " " << j << std::endl;

    std::cout << "Max score: " << max_score << std::endl;
    std::cout << "Alignments:" << std::endl << align1 << std::endl << align2 << std::endl;
    std::cout << "Alignment Length: " << align1.size() << std::endl;

    return align1.size();

}
