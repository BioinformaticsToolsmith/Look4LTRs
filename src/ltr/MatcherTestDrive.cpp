/**
 * The purpose of this program is to test the Matcher class
 */

#include "ScorerTr.h"
#include "Merger.h"
#include "Detector.h"
#include "LtrParameters.h"
#include "Stretch.h"
#include "Element.h"
#include "Node.h"
#include "DirectedGraph.h"
#include "Matcher.h"
#include "Red.h"
#include "../FastaReader.h"

#include <vector>
#include <string.h>
#include <iostream>

int main()
{
        std::string fastaPath = "/home/transposons/Genomes/TAIR10/Fasta/TAIR10.chr1.fa";
        std::string redDir = "/home/transposons/Genomes/TAIR10/RScore/";
        std::string fastaName = fastaPath.substr(fastaPath.find_last_of('/'));
        fastaName = fastaName.substr(0, fastaName.find_last_of('.'));
        FastaReader fr(fastaPath, 1000);
        auto chrom = fr.read()->at(0);
        auto seq = chrom.second;

        // Scoring
        std::cout << "Scoring" << std::endl;
        ScorerTr st(*seq, 13, 250, 34000);

        std::cout << "Forward Merge" << std::endl;

        // Preprocesssing
        Merger forwardM(st.getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                        LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);

        std::cout << "Backward Merge" << std::endl;
        Merger backwardM(st.getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                         LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);

        // forwardM.printStretches();
        // std::cout << std::endl;
        // backwardM.printStretches();
        // std::cout << std::endl;

        auto fStretchVec = forwardM.getStretchVec();
        auto bStretchVec = backwardM.getStretchVec();

        std::cout << "Detecting" << std::endl;
        // Applying Classifier
        Detector dt;
        auto fElementVec = dt.apply(*fStretchVec);
        auto bElementVec = dt.apply(*bStretchVec);

        // dt.printElements(fElement, 10);
        // std::cout << std::endl;
        // dt.printElements(bElement, 10);
        // std::cout << std::endl;

        // Getting Red
        std::cout << "Loading Red" << std::endl;
        Red red;
        red.loadScores(redDir + fastaName + ".scr", chrom);

        std::cout << "Matching" << std::endl;
        // Matching
        Matcher mat{fElementVec, bElementVec, red};
        auto graph = mat.getGraph();
        graph.checkGraph();


        //graph.print("/home/transposons/Projects/Identity/testoutput/nobackward.txt");
        // std::cout << graph << std::endl << std::endl;
        // std::cout << std::endl;

        return 0;
}
