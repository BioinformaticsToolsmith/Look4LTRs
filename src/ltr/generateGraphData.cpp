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
#include <fstream>

int main(int argc, char *argv[]) {

    std::string fastaDir;
    std::string outPath;
    std::string configFile;
    int pa;
    std::string helpMessage = 
        "Example of expected input: generateGraphData -fasta /###/###/Fasta/ -out /Output/ -pa 4\n"
        "--fasta/-f: fasta file directory; train Red on this directory\n"
        "--out/-o  : REQUIRED; output directory; all output files go here\n"
        "--config/-c: configuration file\n"
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

        else if (arg == "--fasta" || arg == "-f") {
            i++;
            fastaDir = argv[i];
        }
        else if (arg == "--out" || arg == "-o") {
            i++;
            outPath = argv[i];
        }
        else if (arg == "--config" || arg == "-c") {
            i++;
            configFile = argv[i];
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

    std::vector<std::string> fastaVec;
    LtrUtility::collectFastaFiles(fastaVec, fastaDir);

    LtrParameters::updateParameters(configFile, false);


    std::unique_ptr<Red> red = std::make_unique<Red>(fastaDir, pa);

        // Database path for Identity training
    std::string dbPath = outPath + "/db.fasta";
    if (std::filesystem::exists(dbPath)) {
        std::filesystem::remove(dbPath);
    }

    // Contains all of the Modules
    std::unordered_map<std::string, ModulePipeline*> moduleMap;

    // Finding Elements
    std::cout << std::endl;
    std::cout << "Finding Repetitive Elements..." << std::endl;
    #pragma omp parallel for schedule(static) num_threads(pa)
    for(int i = 0; i < fastaVec.size(); i++) {
        std::string fastaPath = fastaVec.at(i);
        std::string fastaName = LtrUtility::getFileName(fastaPath);
        std::string chromOut = outPath + "/" + fastaName;


        #pragma omp critical
        {
            std::cout << "Parsing " + fastaName << std::endl;
        }

        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();
        for (auto chrom : *block) {

            #pragma omp critical 
            {
                moduleMap[chromOut] = new ModulePipeline{*red};
            }
            ModulePipeline &mp = *moduleMap[chromOut];

            // Scoring, Merging, and Detecting elements
            mp.buildElements(chrom.second);

        }
        FastaReader::deleteBlock(block);

    }

    std::cout << "Writing Elements to Database..." << std::endl;
    for (int i = 0; i < fastaVec.size(); i++) {
        std::string fastaPath = fastaVec.at(i);
        std::string fastaName = LtrUtility::getFileName(fastaPath);
        std::string chromOut = outPath + "/" + fastaName;


        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();
        for (auto chrom : *block) {

            // Writing elements to database for Identity training
            ModulePipeline &mp = *moduleMap[chromOut];
            mp.writeToDB(dbPath, i, fastaName, chrom.second);
        }
        FastaReader::deleteBlock(block);
    }

    // Training Identity
    std::cout << "Training Identity..." << std::endl;
    auto icStandard = LtrUtility::buildCalculator(LtrParameters::MIN_IDENTITY, dbPath, pa, false);
    auto icRecent = LtrUtility::buildCalculator(LtrParameters::MIN_IDENTITY_RECENT, dbPath, pa, true);
    std::filesystem::remove(dbPath);

    // Looking for LTR RTs
    std::cout << "Building Graphs and writing to file..." << std::endl;
    #pragma omp parallel for schedule(static) num_threads(pa)
    for(int i = 0; i < fastaVec.size(); i++) {
        std::string fastaPath = fastaVec.at(i);
        std::string fastaName = LtrUtility::getFileName(fastaPath);
        std::string chromOut = outPath + "/" + fastaName;

        #pragma omp critical
        {
            std::cout << "Parsing " + fastaPath << std::endl;
        }

        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();

        for (auto chrom : *block) {

            ModulePipeline &mp = *moduleMap[chromOut];


            mp.matchElements(*icStandard, *icRecent, chrom.second);

            // Writing graph to file
            auto graph = mp.getGraph();

            std::string graphPath = chromOut + ".gph";
            std::ofstream graphFile(graphPath);

            graph->write(graphFile);

            graphFile.close();


            delete moduleMap[chromOut];

            break;

        }
        FastaReader::deleteBlock(block);


    }
}