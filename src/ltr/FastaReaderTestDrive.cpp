#include "../FastaReader.h"

#include <iostream>

int main() {

    FastaReader fr("/home/transposons/Projects/Identity/cases/fastareader/test.fa", 10000);
    std::string *x = fr.read()->at(0).first;

    std::cout << *x << std::endl;

    return 0;
}