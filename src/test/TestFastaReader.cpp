#include "../FastaReader.h"

#include <string>
#include <iostream>

int main() {
    std::string fastaPath = "/home/transposons/Genomes/TAIR10/Fasta/TAIR10.chr1.fa";

    FastaReader fr{fastaPath, 100};
    std::cout << *fr.read()->at(0).first << std::endl;
}