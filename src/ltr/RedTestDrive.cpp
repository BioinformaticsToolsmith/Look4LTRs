#include "../red/Red.h"
#include "../FastaReader.h"

#include <vector>
#include <string.h>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <thread>
#include <unordered_map>
#include <sstream>
#include <tuple>
#include <assert.h>

int main(int argc, char *argv[]) {

    std::string fastaDir;
    std::string outputDir;
    if (argc > 1) {
        fastaDir.assign(argv[1]);
        outputDir.assign(argv[2]);
    }
    else {
        fastaDir = "/home/transposons/Genomes/TAIR10/Fasta/";
        outputDir = "/home/transposons/Projects/Identity/testoutput/";
    }

    std::unordered_map<std::string, std::string> seqFileMap;

    int coreCount = std::thread::hardware_concurrency();

    std::cout << "Train RED" << std::endl;
    Red red{fastaDir};    
    
    std::vector<std::filesystem::directory_entry> fastaVec;
    for (auto const& fastaEntry : std::filesystem::directory_iterator{fastaDir}) {
        fastaVec.push_back(fastaEntry);
    }

        #pragma omp parallel for schedule(static) num_threads(coreCount)
    for(int i = 0; i < fastaVec.size(); i++) {
        std::string fastaPath = fastaVec.at(i).path().string();
        std::string fastaName = fastaPath.substr(fastaPath.find_last_of('/') + 1);
        fastaName = fastaName.substr(0, fastaName.find_last_of('.'));    
        std::cout << "Parsing " + fastaPath << std::endl;

        FastaReader fr(fastaPath, 1000);
        for (auto &chrom : *fr.read()) {
            seqFileMap[fastaName + *chrom.first] = outputDir + fastaName + "_" + chrom.first->substr(1);
            auto repeatVec = red.predictRepeats(*chrom.second);

            ofstream file{seqFileMap[fastaName + *chrom.first] + ".rpts"};
            for (auto x : *repeatVec) {
                file << x->getStart() << "\t" << x->getEnd() + 1 << std::endl;
            }

            file.close();
        }
    }
}