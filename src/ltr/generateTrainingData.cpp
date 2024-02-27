#include "ModulePipeline.h"
#include "LtrUtility.h"

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
#include "OutputCpx.h"
#include "Filter.h"
#include "PostProcess.h"

#include "../FastaReader.h"
#include "../IdentityCalculator.h"
#include "../SynDataGenerator.h"
#include "../KmerHistogram.h"
#include "../red/Red.h"
#include "../Matrix.h"
#include "../utility/ILocation.h"

#include <string>
#include <vector>
#include <iostream>
#include <filesystem>
#include <memory>
#include <unordered_set>
#include <algorithm>
#include <map>

std::vector<std::tuple<int, int, int, int>>*  getLTRBedLocations(std::string &bedPath) {

    // left LTR start, left LTR end, right LTR start, right LTR end
    std::vector<std::tuple<int, int, int, int>> *locations = new std::vector<std::tuple<int, int, int, int>>();

    std::ifstream bedFile(bedPath);
    std::string line;
    // Skip the header line
    std::getline(bedFile, line);
    while (std::getline(bedFile, line)) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(line);
        while (std::getline(tokenStream, token, '\t')) {
            tokens.push_back(token);
        }
        int leftStart = std::stoi(tokens.at(3));
        int leftEnd = std::stoi(tokens.at(4));
        int rightStart = std::stoi(tokens.at(5));
        int rightEnd = std::stoi(tokens.at(6));

        locations->push_back(std::make_tuple(leftStart, leftEnd, rightStart, rightEnd));
    }
    return locations;
    // up to user to delete
}

int main(int argc, char *argv[]) {

    std::string fastaRedPath;
    std::string syntheticPath;
    std::string bedPath;
    std::string outPath;
    int pa;
    std::string helpMessage = 
        "Example of expected input: generateTrainingdata --synthetic /###/###/SyntheticFasta/ -fasta /###/###/Fasta/ --bed /###/###/Bed/ --out /Output/ -pa 4\n"
        "--synthetic/-s: synthetic fasta file directory; generate stretches from these files\n"
        "--fasta/-f: fasta file directory; train Red on this directory\n"
        "--bed/-b: bed file directory; get repeat content from interiors of LTR-retrotransposons using these files\n"
        "--out/-o  : output directory; all output files go here\n"
        "--parallel/-pa   : number of cores to use\n"
        "--help/-h : prints out this message\n\n";

    if (argc < 2) {
        std::cout << helpMessage;
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            std::cout << helpMessage;
            return 0;
        }
        else if (arg == "--synthetic" || arg == "-s") {
            i++;
            syntheticPath = argv[i];
        }
        else if (arg == "--fasta" || arg == "-f") {
            i++;
            fastaRedPath = argv[i];
        }
        else if (arg == "--bed" || arg == "-b") {
            i++;
            bedPath = argv[i];
        }
        else if (arg == "--out" || arg == "-o") {
            i++;
            outPath = argv[i];
        }
        else if (arg == "--parallel" || arg == "-pa") {
            i++;
            pa = std::stoi(argv[i]);
        }
        else {
            std::cout << "Invalid argument: " << arg << std::endl;
            std::cout << helpMessage;
            return 1;
        }
    }

    std::vector<std::string> syntheticVec;
    LtrUtility::collectFastaFiles(syntheticVec, syntheticPath);

    std::vector<std::string> fastaVec;
    LtrUtility::collectFastaFiles(fastaVec, fastaRedPath);

    if (syntheticVec.size() != fastaVec.size()) {
        std::cout << "The number of synthetic fasta files and fasta files must be the same" << std::endl;
        return 1;
    }

    std::vector<std::string> bedVec;
    LtrUtility::collectBedFiles(bedVec, bedPath);


    if (bedVec.size() != fastaVec.size()) {
        std::cout << "The number of bed files and fasta files must be the same" << std::endl;
        return 1;
    }

    std::sort(syntheticVec.begin(), syntheticVec.end(), [](const std::string& a, const std::string& b) { return LtrUtility::getFileName(a) < LtrUtility::getFileName(b); });

    std::sort(fastaVec.begin(), fastaVec.end(), [](const std::string& a, const std::string& b) { return LtrUtility::getFileName(a) < LtrUtility::getFileName(b); });

    std::sort(bedVec.begin(), bedVec.end(), [](const std::string& a, const std::string& b) { return LtrUtility::getFileName(a) < LtrUtility::getFileName(b); });

    std::unique_ptr<Red> red = std::make_unique<Red>(fastaRedPath, pa);


    /*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*
     * Getting Training Data for Detector Module
    *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*/
    for (int i = 0; i < syntheticVec.size(); i++) {

        std::string syntheticFastaPath = syntheticVec.at(i);
        std::string syntheticFastaName = LtrUtility::getFileName(syntheticFastaPath);
        std::string semiChromOut = outPath + "/" + syntheticFastaName;


        std::string fastaPath = fastaVec.at(i);
        std::string fastaName = LtrUtility::getFileName(fastaPath);

        #pragma omp critical
        {
            std::cout << "Parsing " + fastaName << std::endl;
        }

        FastaReader sfr(syntheticFastaPath, 1000);
        FastaReader fr(fastaPath, 1000); 
        auto semiBlock = sfr.read();
        auto block = fr.read();

        for (size_t i = 0; i < block->size(); ++i) {
            // Access elements using index i
            auto& semiChrom = (*semiBlock)[i];
            auto& chrom = (*block)[i];

            ModulePipeline *mpPtr = new ModulePipeline{*red};

            // Scoring and merging stretches
            mpPtr->buildStretches(semiChrom.second);

            // Writing feature matrix to file
            mpPtr->writeFeatureMatrix(semiChromOut + ".fmx", semiChrom.second, chrom.second);

            // writing stretches to file
            mpPtr->writeStretches(semiChromOut + ".fs", semiChromOut + ".bs");

            delete mpPtr;

            break;
        }
        FastaReader::deleteBlock(semiBlock);
        FastaReader::deleteBlock(block);
    }

    /*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*
     * Getting Training Data for Interior repetition check
    *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*/

    for (int i = 0; i < fastaVec.size(); i++) {

        std::string fastaPath = fastaVec.at(i);
        std::string fastaName = LtrUtility::getFileName(fastaPath);
        std::string chromOut = outPath + "/" + fastaName;


        std::string bedPath = bedVec.at(i);
        std::string bedName = LtrUtility::getFileName(bedPath);

        std::vector<std::tuple<int, int, int, int>> *ltrLocations = getLTRBedLocations(bedPath);

        #pragma omp critical
        {
            std::cout << "Parsing " + fastaName << std::endl;
        }

        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();

        // open a file at chromOut + ".rsc" and write red scores to it
        std::ofstream redScoreFile(chromOut + ".rsc");

        for (size_t i = 0; i < block->size(); ++i) {
            // Access elements using index i
            auto& chrom = (*block)[i];

            std::string *seq = chrom.second;

            for (auto &loc : *ltrLocations) {
                int leftEnd = std::get<1>(loc);
                int rightStart = std::get<2>(loc);

                int interiorSize = rightStart - leftEnd;

                std::string interior = seq->substr(leftEnd, interiorSize);

                auto scoreVec = red->score(interior);

                // on one line, write the scores for the interior down. 
                for (size_t i = 0; i < scoreVec.size(); ++i) {
                    redScoreFile << scoreVec[i];
                    if (i != scoreVec.size() - 1) {
                        redScoreFile << " ";
                    }
                }
                redScoreFile << std::endl;
            }

            break;
        }
        delete ltrLocations;
        redScoreFile.close();
        FastaReader::deleteBlock(block);
    }


}
