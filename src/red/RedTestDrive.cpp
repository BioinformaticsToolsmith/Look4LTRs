#include "Red.h"

#include <string>
#include <iostream>

int main() {
    /**
     * Parameters as follows:
     * genome directory 
     * core count 
     * k 
     * order 
     * gaussian mask half width
     * threshold
    */
    std::string goodPath = "/home/transposons/Genomes/TAIR10/Fasta/";
    std::string badPath = "/home/transposons/Genomes/TAIR10/GNM"; // DOES NOT EXIST 

    int goodCore = 4;
    int badCore = -1;

    int goodK = 13;
    int badK = -15;
    
    int goodOrder = 2;
    int badOrder = -3;
    
    double goodGau = 30.0;
    double badGau = -40.0;

    double goodThr = 3.0;
    double badThr = -1.0;

    int goodMin = 3;
    int badMin = -4;

    std::string goodCnd = "/home/transposons/Projects/Identity/testoutput/";
    std::string badCnd = "/home/transposons/Projects/Identity/testoutput/doesntexist/";

    std::string goodTblPath = "/home/transposons/Projects/Identity/testoutput/table.txt";
    std::string badTblPath = "/home/transposons/Projects/Identity/doesntexist/table.txt";

    std::string goodHmmPath = "/home/transposons/Projects/Identity/testoutput/hmm.txt";
    std::string badHmmPath = "/home/transposons/Projects/Identity/doesntexist/hmm.txt";

    std::string goodRptPath = "/home/transposons/Projects/Identity/testoutput/";
    std::string badRptPath = "/home/transposons/Projects/Identity/doesntexist/";

    int goodFormat = 2;
    int badFormat = -5;

    std::string goodScoPath = "/home/transposons/Projects/Identity/testoutput/";
    std::string badScoPath = "/home/transposons/Projects/Identity/doesntexist/";

    std::string goodMskPath = "/home/transposons/Projects/Identity/testoutput/";
    std::string badMskPath = "/home/transposons/Projects/Identity/doesntexist/";

    std::string goodDirPath = "/home/transposons/Genomes/Test/";
    std::string badDirPath = "/home/transposons/Genomes/DoesntExist!/";

    Red red{goodPath, {}, {}, {}, {}, {}, {}, {}};
    // red.printTable(badTblPath);
    // red.printModel(badHmmPath);

    //red.scan(goodRptPath, goodFormat, goodScoPath, goodMskPath, goodDirPath);

    // This is inside of Arabidopsis thaliana
    std::string ATCopia58LTR = "TGTTGAAAGTTAAACTTGATTTTGAATCAAGTTTAATTATTGGATCAATTATCCAATAATTAATTATGGCCAAATCCAAGTTCTAGAGTTTTCTCTAGAAATATCATCATTTCCACCTCCTTAAAAGATTCTAGAAATTTTCTAGAATCATCTTCCACCTCCTTAAACATAAAAATCTAGATACTCTAATAGAATAATCTAGATAATTTGAATAATGTAATCTAGATCTTATGTAAGAACTCTCTAGACTTAGGATTAAAATATTTTAGATATTTTGTAGTTTGGAGGCTATAAATACCTCCTCCCCCTCTCAAATGTTGCAATGTTGTGAAGTTGTATTCAAGTTTAAAGCAAAGTAATAAAAGTTCTATTTCCTAAAAAACTCTCTCAAAACACTTAAACACTTTCTCCATTACCTCTAAAAGAATTTTACTCTAACA";
    std::string nonRepetitive = "AAGATCTTTCATCATCCTCCAGGAAAAGAACGTTTTAAAATTTTATATTCCAGAAGAAAACAAACACTTTTATATTGTGTCGTTGAGGTTGAGTTGTGTTTGGAAGATAAGTTTATTGACCTATTGATCTGTAACTTCATAAGATTTTGAAACGTTAGAAGATTCAAAAGAATGTTTGTTTTGCATTTTTTTTAAAAAATATTGCAAAAGAATATATAGATTCTATTTTAGATAACTATCATCAAATAGTTTTAAAGAAAAAAACAAAAATTCTGTTTCTTGAATGAATGACTTTTAGATCTTATTTTTCGCCTTTTGCAGGATTGAAAAATACTAAAAGATAAAAAGATTTGAATCTTGCACTTTAGCCCCAGAATTACCGTTTTTGATGATTTTGATGCAAAAGTGAGACTACTAAACTTTTGAGAGT";

    auto scoreVec = red.score(nonRepetitive);
    int count = 0;
    for (auto x : scoreVec) {
        if (x != 0) {
            count += 1;
        }
    }

    std::cout << static_cast<double>(count) / scoreVec.size() << std::endl; 

    

}