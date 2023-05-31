#include "ScorerTr.h"
#include "Merger.h"
#include "Detector.h"
#include "LtrParameters.h"
#include "Stretch.h"
#include "Element.h"
#include "DirectedGraph.h"
#include "Matcher.h"
#include "OutputBed.h"
#include "OutputRtr.h"
#include "Filter.h"
#include "PostProcess.h"
#include "../Matrix.h"


#include "../FastaReader.h"
#include "../IdentityCalculator.h"
#include "../SynDataGenerator.h"
#include "../KmerHistogram.h"
#include "../red/Red.h"

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

void writeStretches(std::string filePath, std::vector<Stretch> &stretchVec);


int main(int argc, char *argv[]) {

    std::string fastaDir(argv[1]);
    std::string outputDir(argv[2]);
    std::unordered_map<std::string, std::string> seqFileMap;
    int coreCount = std::thread::hardware_concurrency();

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
            int id = 1;
            seqFileMap[fastaName + *chrom.first] = outputDir + fastaName + "_" + chrom.first->substr(1);

            // Scoring
            std::cout << "Scoring" << *chrom.first << std::endl;
            ScorerTr st(*chrom.second, 13, 250, 34000);
            // st.printForwardScores(seqFileMap[fastaName + *chrom.first] + ".fscr");
            // st.printBackwardScores(seqFileMap[fastaName + *chrom.first] + ".bscr");

            // Forward Merging
            std::cout << "Forward Merge" << *chrom.first << std::endl;
            Merger forwardM(st.getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                            LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);

            // Backward Merging
            std::cout << "Backward Merge" << *chrom.first <<  std::endl;
            Merger backwardM(st.getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                            LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);
            

            std::cout << "Creating stretch files..." << std::endl;
            writeStretches(outputDir + "/" + fastaName + "_Forward", *forwardM.getStretchVec());
            writeStretches(outputDir + "/" + fastaName + "_Backward", *backwardM.getStretchVec());

        }
    }

    return 0;
  
}

void writeStretches(std::string filePath, std::vector<Stretch> &stretchVec) {

    Matrix stretchMatrix{stretchVec.size(), 3};
    for (int i = 0; i < stretchVec.size(); i++) {
        stretchMatrix(i, 0) = stretchVec.at(i).getStart();
        stretchMatrix(i, 1) = stretchVec.at(i).getEnd();
        stretchMatrix(i, 2) = stretchVec.at(i).getMedianHeight();
    }

    stretchMatrix.printToFile(filePath + ".stc");

}