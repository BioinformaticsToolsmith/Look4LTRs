#include "../IdentityCalculator.h"
#include "../SynDataGenerator.h"

#include <string>
#include <thread>
#include <iostream>

template <typename T>
void start(SynDataGenerator *, int, double, bool, bool);
int main()
{

    std::string dbPath = "/home/transposons/Data/db.fasta";
    int coreCount = std::thread::hardware_concurrency();
    SynDataGenerator dg{dbPath, 0.8, coreCount};

    int64_t maxLength = dg.getMaxLength();

    // Determine histogram data type
    if (maxLength <= std::numeric_limits<int8_t>::max())
    {
        start<int8_t>(&dg, coreCount, 0.8, true, true);
    }
    else if (maxLength <= std::numeric_limits<int16_t>::max())
    {
        start<int16_t>(&dg, coreCount, 0.8, true, true);
    }
    else if (maxLength <= std::numeric_limits<int32_t>::max())
    {
        start<int32_t>(&dg, coreCount, 0.8, true, true);
    }
    else if (maxLength <= std::numeric_limits<int64_t>::max())
    {
        start<int64_t>(&dg, coreCount, 0.8, true, true);
    }
    else
    {
        std::cout << "ReaderAlignerCoordinator warning: ";
        std::cout << "Overflow is possible however unlikely.";
        std::cout << std::endl;
        std::cout << "A histogram entry is 64 bits." << std::endl;
        start<int64_t>(&dg, coreCount, 0.8, true, true);
    }
}

template <typename T>
void start(SynDataGenerator *dg, int core, double thresh, bool skip, bool relax)
{
    KmerHistogram<uint64_t, T> kTable(dg->getK());
    KmerHistogram<uint64_t, uint64_t> monoTable(1);

    std::string one = "GTCGTCGAACAACAATTTCATGCCACCCGTCCTCCACCCATGTTTTGGCCACTAAGACAGGTAAGGCTGTTTTTGGCCTCGCGTGAGCTACAACACATGTTTTCATGGCCGAACAACCAATTTAGTGTCCAACCATAGTACACTAGTGTTGAAACCATAGTACACATTTTTGTCCCTGGAGGCCTGTAAGGCTATTTTTGGCCTCCTGCGACACATGTTTTCTTCGTCGAACAACAATTTCATGCCTCCCGCCCGCCAAACATGTTTTGGCTATTGACATGGGTAAGGCTGTTTATAGCATCTGTTGAGCTACATCACACAAAACACCTAAATCCTAAACCCCGAGCCCCAAACCCTAAACCATGAACCCGGAACCGCGAACCCTTAGACTAAAACCCGACCCCCAAAACACAAAATCACAAATCCCAAATTCCAAACCCTAAACATTAAACCCCAAACCCTAACCCTCAAATCGAAACCTGAAATCTCGAACCCCAAAAGCCTAAGCTCTAAACCCCGAGCCCCAAACCCAAAAGCCTGAACATTGAACCCCGAACCCTAATTTGTGAACCCTATACCCCGAACCCTGAAACAAAGACTCGACCACCAAAACACAAAACTCCAAACCCCAAACTCTAAACCCTAAAACACAACCCTAAAAAGTAAAGTAAGGTTGTTTTTGGCCTCGCGTGAGCTACAACACACGTTTTCATGGCCGAACAACCAATTTAGTGTCCAAACCATAGTACACTAGTGTTCGAACCATTGTACACGTATTTGGCCCTGGAGGCCGGTAAGGCTATTTTT";
    // std::string two = "AGGTGTACGAGCCTCTGGTCGATGATCAATGGCAACACAACCCCATTTTTGTCGATAATAGCCATGAATGATCATTTTCAATAATACCGAAGGCTAACACCTATAGATT";
    std::string two = one;

    IdentityCalculator<T> ic{dg, core, thresh, skip, relax};
    ic.score(&one, &two);

    // T *kHistOne = kTable.build(&one);
    // uint64_t *monoHistOne = monoTable.build(&one);
    // T *kHistTwo = kTable.build(&two);
    // uint64_t *monoHistTwo = monoTable.build(&two);

    // int oneLen = one.size();
    // int twoLen = two.size();

    // double score = ic.score(kHistOne, kHistTwo, monoHistOne, monoHistTwo, ic.calcRatio(oneLen, twoLen), oneLen, twoLen);
    // std::cout << score << std::endl;

    // delete[] kHistOne;
    // delete[] monoHistOne;
    // delete[] kHistTwo;
    // delete[] monoHistTwo;
}

double calcIdentity(std::string seq1, std::string seq2)
{
    double id = 0.0;

    return id;
}