#include "OutputBed.h"
#include "OutputRtr.h"
#include "../FastaReader.h"
#include "RTSolo.h"
#include "RTComplete.h"
#include "RT.h"
#include "Element.h"
#include "Stretch.h"

#include <vector>
#include <iostream>
#include <string>

RT* makeRTComplete(int l_start, int l_end, int r_start, int r_end);
RT* makeRTSolo(int start, int end);

int main() {
    std::string outputDir = "/home/transposons/Projects/Identity/testoutput/";
    std::string fastaPath = "/home/transposons/Projects/Identity/cases/fastareader/test.fa";

    OutputBed oBed{outputDir, fastaPath};
    OutputRtr oRtr{outputDir, fastaPath};
    FastaReader fr(fastaPath, 10000);
    for (auto &chrom : *fr.read()) {

        std::vector<RT*> rtVec{makeRTComplete(0, 100, 900, 1000), 
                                makeRTSolo(1050, 1150), 
                                makeRTComplete(1200, 1300, 2600, 2700), 
                                makeRTComplete(1500, 1800, 2200, 2250)};

        rtVec[2]->nest(rtVec[3], true);

        oBed.write(*chrom.first, rtVec);
        oRtr.write(*chrom.first, rtVec);

        for (auto ptr : rtVec) {
            delete ptr;
        }
    }

    return 0;
}


RT* makeRTComplete(int l_start, int l_end, int r_start, int r_end) {
    RT* r;

    Stretch* left = new Stretch{l_start, l_end, Stretch::K, true};
    Stretch* right = new Stretch{r_start, r_end, Stretch::K, false};

    r = new RTComplete{new Element{*left}, new Element{*right}};

    return r;
}

RT* makeRTSolo(int start, int end) {
    RT* r;

    Stretch* l = new Stretch{start, end, Stretch::K, true};

    r = new RTSolo{new Element{*l}};

    return r;
}


