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
#include "LtrUtility.h"

#include "ModulePipeline.h"

#include "../FastaReader.h"
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

int main(int argc, char*argv[]) {

    std::vector<std::string> args(argv + 1, argv + argc);

    std::vector<std::string> fastaArgVec;
    std::vector<std::string> trainArgVec;
    std::string outPath = "";
    std::string configPath = "";
    int pa = 1;
    bool help = false;
    std::string helpMessage = 
        "Example of expected input: look4ltrs -fasta /###/###/Fasta/ /###/####/Fasta -out /Output/ -pa 4\n"
        "--fasta/-f: fasta file directories; variable number of arguments; trains and predict on these genomes\n"
        "--train/-t: fasta file directories; variable number of arguments; trains on these genomes, no prediction\n"
        "--out/-o  : REQUIRED; output directory; all output files go here\n"
        "--config/-c : library stats for detector module (optional)\n"
        "--parallel/-pa   : number of cores to use\n"
        "--help/-h : prints out this message";


    /**
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     *                  ARGUMENT PARSING
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     * Refer to the help message above for what these arguments mean
    */    
    int i = 0;
    while (i < args.size()) {
        if (args[i] == "--fasta" || args[i] == "-f") {
            i++;
            while (i < args.size() && args[i][0] != '-') {
                fastaArgVec.push_back(args[i] + "/");
                i++;
            }
        }
        else if (args[i] == "--train" || args[i] == "-t") {
            i++;
            while (i < args.size() && args[i][0] != '-') {
                trainArgVec.push_back(args[i] + "/");
                i++;
            }
        }
        else if (args[i] == "--out" || args[i] == "-o") {
            i++;
            outPath = args[i] + "/";
            i++;
        }
        else if (args[i] == "--config" || args[i] == "-c") {
            i++;
            configPath = args[i];
            i++;
        }
        else if (args[i] == "--parallel" || args[i] == "-pa") {
            i++;
            pa = std::stoi(args[i]);
            pa = pa == -1 ? std::thread::hardware_concurrency() : pa;
            i++;
        }
        else if (args[i] == "--help" || args[i] == "-h") {
            help = true;
            i++;
        }
        else {
            std::cerr << "Invalid argument: " << args[i] << std::endl;
            return 1;
        }
    }


    /**
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     *                  ARGUMENT VALIDATION
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    */    

    // Regardless if other arguments are passed, if help is given then print the help message and return
    // If the user only provided the name of the executable, do the same
    if (help || argc == 1) {
        std::cout << helpMessage << std::endl;
        return 0;
    }

    // The output directory is required; There is far too much information to print to terminal.
    if (outPath.empty()) {
        std::cerr << "--out is a necessary parameter!" << std::endl;
        return 1;
    }

    // We need something to train and predict on, after all
    if (fastaArgVec.empty()) {
        std::cerr << "Must pass at least one fasta directory to --fasta!" << std::endl;
        return 1;
    }


    // Checking fasta directories
    for (auto& fastaDir : fastaArgVec) {
        if (!std::filesystem::exists(fastaDir)) {
            std::cerr << fastaDir << " does not exist!" << std::endl;
            std::cerr << "Error in argument --fasta" << std::endl;
            return 1;
        }
    }

    // Checking train directories
    for (auto& fastaDir : trainArgVec) {
        if (!std::filesystem::exists(fastaDir)) {
            std::cerr << fastaDir << " does not exist!" << std::endl;
            std::cerr << "Error in argument --train" << std::endl;
            return 1;
        }
    }


    // Checking output directory
    if (!std::filesystem::exists(outPath)) {
        std::cerr << outPath << " does not exist!" << std::endl;
        std::cerr << "Error in argument --out" << std::endl;
        return 1;
    }

    // Checking core value
    if (pa < 1) {
        std::cerr << "Invalid number of cores: " << pa << std::endl;
        return 1;
    }
    if (pa > std::thread::hardware_concurrency()) {
        std::cerr << "Number of cores exceeds hardware concurrency!" << std::endl;
        return 1;
    }

    /**
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     *            COLLECTING FASTA FILES
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    */

    // Collecting training fasta files;
    // Used only for training RED and Identity
    std::vector<std::string> trainVec;
    LtrUtility::collectFastaFiles(trainVec, trainArgVec);

    // Collecting prediction fasta files;
    // Trains RED, Identity, and is predicted upon
    std::vector<std::string> predVec;
    LtrUtility::collectFastaFiles(predVec, fastaArgVec);

    
    // Making sure there aren't duplicates in the passed in files. That would make our fasta file tracking system go out of wack. 

    // Checking duplicates in the training files
    std::unordered_set<std::string> trainSet(trainVec.begin(), trainVec.end());
    if (trainVec.size() != trainSet.size()) {
        std::cerr << "Duplicate fasta files were passed into  --train! Check your file names!" << std::endl;
        return 1;
    }

    // Checking duplicates in the predicting files
    std::unordered_set<std::string> predSet(predVec.begin(), predVec.end());
    if (predVec.size() != predSet.size()) {
        std::cerr << "Duplicate fasta files were passed into --fasta! Check your file names!" << std::endl;
        return 1;
    }

    // Retrieving all fasta files from the training and predicting sets
    std::unordered_set<std::string> fastaSet(trainVec.begin(), trainVec.end());
    fastaSet.insert(predVec.begin(), predVec.end());
    std::vector<std::string> fastaVec(fastaSet.begin(), fastaSet.end());
    std::sort(fastaVec.begin(), fastaVec.end());



    // Ensuring user sees what files are being used
    std::cout << "Training-only files:" << std::endl;
    for (auto& fastaFile : trainVec) {
        std::cout << fastaFile << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Training-and-Prediction files:" << std::endl;
    for (auto& fastaFile : predVec) {
        std::cout << fastaFile << std::endl;
    }
    std::cout << std::endl;

    /**
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     *                  MAIN CODE
     * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    */

    // Updating detector vectors / library statistics if provided
    if (!configPath.empty()) {
        LtrParameters::updateParameters(configPath);
    }

    for (auto& val : LtrParameters::MEAN_VECTOR) {
        std::cout << val << std::endl;
    }

    // Training RED
    std::string redTrainPath;
    std::vector<std::filesystem::path> symlinkVec;
    std::unique_ptr<Red> red;

    // Do we need to build symbolic links? Yes if there is more than one directory to consider
    bool needSymlinks = fastaArgVec.size() > 1 || trainArgVec.size() > 0;

    // Building symbolic links if RED has to train on files from different directories.
    try {
        // Different directories?
        if (needSymlinks) {
            redTrainPath = outPath + "/temp/";
            if (std::filesystem::exists(redTrainPath)) {
                std::filesystem::remove_all(redTrainPath);
            }
            // Building symbolic links
            std::filesystem::create_directory(redTrainPath);
            for (auto& fastaFile : fastaVec) {
                std::filesystem::path file(fastaFile);
                std::filesystem::path symlinkPath = redTrainPath / file.filename();
                std::filesystem::create_symlink(file, symlinkPath);
                symlinkVec.push_back(symlinkPath);
            }
        }
        // Just one directory, no need for messy, symbolic-link management
        else {
            redTrainPath = fastaArgVec.at(0);
        }

        // Training Red right here
        red = std::make_unique<Red>(redTrainPath, pa);
    }
    // If we fail somewhere above, remove any created symbolic link.
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        for (auto& symlink : symlinkVec) {
            std::filesystem::remove(symlink);
        }
        return 1;
    }

    // After training RED, we don't need the symbolic links
    for (auto& symlink : symlinkVec) {
        std::filesystem::remove(symlink);
    }
    // If we created the temporary directory for symbolic links, remove it
    if (needSymlinks) {
        std::filesystem::remove_all(redTrainPath);
    }

    // Database path for Identity training
    std::string dbPath = outPath + "db.fasta";
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

        #pragma omp critical
        {
            std::cout << "Parsing " + fastaPath << std::endl;
        }

        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();
        for (auto chrom : *block) {

            std::string chromName = chrom.first->substr(1);
            std::string chromOut = outPath + fastaName + "_" + chromName;

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

        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();
        for (auto chrom : *block) {
            std::string chromName = chrom.first->substr(1);
            std::string chromOut = outPath + fastaName + "_" + chromName;

            // Writing elements to database for Identity training
            ModulePipeline &mp = *moduleMap[chromOut];
            mp.writeToDB(dbPath, i, chromName, chrom.second);
        }
        FastaReader::deleteBlock(block);
    }



    // Training Identity
    std::cout << "Training Identity..." << std::endl;
    auto icStandard = LtrUtility::buildCalculator(LtrParameters::MIN_IDENTITY, dbPath, pa, false);
    auto icRecent = LtrUtility::buildCalculator(LtrParameters::MIN_IDENTITY_RECENT, dbPath, pa, true);
    std::filesystem::remove(dbPath);

    // Looking for LTR RTs
    std::cout << "Looking for LTR RTs..." << std::endl;
    #pragma omp parallel for schedule(static) num_threads(pa)
    for(int i = 0; i < predVec.size(); i++) {
        std::string fastaPath = predVec.at(i);
        std::string fastaName = LtrUtility::getFileName(fastaPath);

        #pragma omp critical
        {
            std::cout << "Parsing " + fastaPath << std::endl;
        }

        // Output file handlers
        OutputBed oBed{outPath + "/Bed", fastaPath};
        OutputRtr oRtr{outPath + "/Rtr", fastaPath};
        OutputCpx oCpx{outPath + "/Cpx", fastaPath};

        FastaReader fr(fastaPath, 1000);
        auto block = fr.read();

        for (auto chrom : *block) {

            std::string chromName = chrom.first->substr(1);
            std::string chromOut = outPath + fastaName + "_" + chromName;

            ModulePipeline &mp = *moduleMap[chromOut];

            mp.matchElements(*icStandard, *icRecent, chrom.second);
            mp.findRTs();

            // Found LTR RTs; overall LTR RT vec
            auto rtVecPtr = mp.getRtVec();
            // The complex LTR RT regions
            auto complexVecPtr = mp.getComplexVec();

            // #pragma omp critical 
            // {
            //     std::cout << "Deep Nests " << chromName << std::endl;
            // }

            // Finding Recently Nested LTR RTs
            mp.findDeepNests(*icStandard, *icRecent, chrom.second);

            // std::cout << "Processing" << std::endl;
            // Re-nest the elements and extend the ends of the two LTRs of each RT by k-1; the scoring module loses out on the last k-1 nucleotides by algorithmic design (big word)
            #pragma omp critical 
            {
                std::cout << "Processing " << chromName << std::endl;
            }
            // Processing the elements, nesting, extending
            mp.process(*icStandard, chrom.second);

            // #pragma omp critical 
            // {
            //     std::cout << "Filter " << chromName << std::endl;
            // }

            // Filtering out the LTR RTs that do not meet the structural features of an LTR RT or exhibit high similarity to other repetitive elements
            mp.filter(*icStandard, chrom.second);


             
            // // Writing to output
            std::cout << "Writing to Output " << *chrom.first << std::endl;
            std::string chromID = chromName.substr(0, chromName.find(" "));
            oBed.write(chromID, *rtVecPtr);
            oRtr.write(chromID, *rtVecPtr);
            oCpx.write(chromID, *complexVecPtr);

            for (auto& r : *rtVecPtr) {
                delete r;
            }
            for (auto& c : *complexVecPtr) {
                delete c;
            }
            delete moduleMap[chromOut];
            moduleMap[chromOut] = nullptr;
        }
        FastaReader::deleteBlock(block);
    }

    for (auto& [key, value] : moduleMap) {
        //value->removeElements();
        if (value != nullptr) {
            delete value;
        }
    }

    std::cout << "Output written to " << outPath << std::endl;
    return 0;
}