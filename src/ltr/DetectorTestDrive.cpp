#include "Detector.h"
#include "Stretch.h"
#include "Element.h"
#include "LtrParameters.h"
#include "Merger.h"
#include "ScorerTr.h"
#include "../FastaReader.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

void writeFile(std::string path, std::vector<Element> eleVec);
int main(int argc, char *argv[]) {

	// std::string input(argv[1]);
	// std::string outDir(argv[2]);
	// std::string fastaName = filesystem::path(argv[1]).stem().string();

    std::string input = "/home/transposons/Projects/Identity/cases/nestedparser/testcasefasta/3elementnest.fa";
    std::string outDir = "/home/transposons/Projects/Identity/cases/nestedparser/testcasefasta/";
    std::string fastaName = "3elementnest";


	FastaReader fr(input, 1000);
	std::string *x = fr.read()->at(0).second;

	ScorerTr st(*x, 13, 250, 34000);

    Merger mf(st.getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
        LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);
    Merger mb(st.getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
        LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);

    // Make Detector object; 
    Detector detector;

    auto fElement = detector.apply(*mf.getStretchVec());
    auto bElement = detector.apply(*mb.getStretchVec());

    writeFile(outDir + "/" + fastaName + ".fele", fElement);
    writeFile(outDir + "/" + fastaName + ".bele", bElement);

}

void writeFile(std::string path, std::vector<Element> eleVec) {
    std::ofstream file;
    file.open(path);
    for (auto &x : eleVec) {
        file << x << std::endl;
    }
    file.close();
}


    // // All but the first pair should merge
    // std::vector<Stretch> a {
    //     makeStretch(10, 200, 400),
    //     makeStretch(800, 890, 1000),
    //     makeStretch(900, 940, 1000),
    //     makeStretch(950, 1000, 978)
    // };

    // std::vector<Stretch> b;
    // std::vector<Stretch> c{makeStretch(10, 200, 500)};
    // std::vector<Stretch> d{makeStretch(10, 200, 500), makeStretch(260, 300, 500)};
    // std::vector<Stretch> e{makeStretch(10, 200, 500), makeStretch(8000, 9000, 1200)};
    // std::vector<Stretch> f{makeStretch(10, 200, 500), makeStretch(210, 300, 510)};
    // std::vector<Stretch> g{makeStretch(10, 200, 500), makeStretch(210, 300, 510), makeStretch(8000, 9000, 1200)};
    // std::vector<Stretch> h{makeStretch(10, 200, 500), makeStretch(210, 300, 510), makeStretch(320, 360, 500)};
