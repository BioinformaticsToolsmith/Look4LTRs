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


#include "../FastaReader.h"
#include "../IdentityCalculator.h"
#include "../SynDataGenerator.h"
#include "../KmerHistogram.h"
#include "../red/Red.h"
#include "../Matrix.h"

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

int writeToFiles(std::string dbPath, std::string elePath, std::vector<Element> &eleVec, std::string *seq, int fastaId, int id);
std::vector<Element> readEleFile(std::string elePath);
void writeEleFile(std::string elePath, std::vector<Element> &eleVec);
void writePrediction(std::string path, std::vector<int> pred);
void writeFeatures(std::string path, Matrix features);
void writeScores(std::string path, std::vector<int> scoreVec, std::string chrom);
void writeStretches(std::string filePath, std::vector<Stretch> &stretchVec);
void writeRepeats(std::string filePath, std::vector<ILocation*>* repeatVec, std::string chrom);

int main(int argc, char *argv[]) {

    std::string fastaDir;
    std::string trainFastaDir;
    std::string outputDir;

    if (argc > 1) {
        fastaDir.assign(argv[1]);
        trainFastaDir.assign(argv[2]);
        outputDir.assign(argv[3]);
    }
    else {
        fastaDir = "/home/transposons/Projects/Identity/Results/NoFilter1/TAIR10/MissingFasta/";
        trainFastaDir = "/home/transposons/Genomes/TAIR10/Fasta/";
        outputDir = "/home/transposons/Projects/Identity/testoutput/";
    }

    std::string dbPath = outputDir + "db.fasta";
    std::unordered_map<std::string, std::string> seqFileMap;
    int coreCount = std::thread::hardware_concurrency();
    std::filesystem::remove(dbPath);

    std::cout << "Generating Red files" << std::endl;
    Red red{trainFastaDir};

    std::vector<std::filesystem::directory_entry> fastaVec;
    for (auto const& fastaEntry : std::filesystem::directory_iterator{fastaDir}) {
        fastaVec.push_back(fastaEntry);
    }
    std::vector<std::filesystem::directory_entry> trainFastaVec;
    for (auto const& fastaEntry : std::filesystem::directory_iterator{trainFastaDir}) {
        trainFastaVec.push_back(fastaEntry);
    }    
    
    int maxLen = 0;
    #pragma omp parallel for schedule(static) num_threads(coreCount)
    for(int i = 0; i < trainFastaVec.size(); i++) {
        std::string fastaPath = trainFastaVec.at(i).path().string();
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

            // Forward Merging
            std::cout << "Forward Merge" << *chrom.first << std::endl;
            Merger forwardM(st.getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                            LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);

            // Backward Merging
            std::cout << "Backward Merge" << *chrom.first <<  std::endl;
            Merger backwardM(st.getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                            LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);
            

            // Applying Classifier
            std::cout << "Detecting" << *chrom.first << std::endl;
            Detector dt{red, *chrom.second};

            auto fElement = dt.apply(*forwardM.getStretchVec());
            auto bElement = dt.apply(*backwardM.getStretchVec());

            for (auto e : fElement) {
                int v = e.getSize();
                maxLen = v > maxLen ? v : maxLen;
            }
            for (auto e : fElement) {
                int v = e.getSize();
                maxLen = v > maxLen ? v : maxLen;
            }

            id = writeToFiles(dbPath, seqFileMap[fastaName + *chrom.first] + ".fele", fElement, chrom.second, i, id);
            id = writeToFiles(dbPath, seqFileMap[fastaName + *chrom.first] + ".bele", bElement, chrom.second, i, id);

            std::filesystem::remove(seqFileMap[fastaName + *chrom.first] + ".fele");
            std::filesystem::remove(seqFileMap[fastaName + *chrom.first] + ".bele");

        }
    }

    std::cout << maxLen << std::endl;
    SynDataGenerator dg{dbPath, LtrParameters::MIN_IDENTITY, coreCount};
    std::filesystem::remove(dbPath);
    int64_t maxLength = dg.getMaxLength();
    assert (maxLength <= std::numeric_limits<int32_t>::max());
    double threshold = 0.8;

    IdentityCalculator<int32_t>* ic = new IdentityCalculator<int32_t>{&dg, coreCount, threshold, true, true};

    //#pragma omp parallel for schedule(static) num_threads(coreCount)
    for(int i = 0; i < fastaVec.size(); i++) {
        std::string fastaPath = fastaVec.at(i).path().string();
        std::string fastaName = fastaPath.substr(fastaPath.find_last_of('/') + 1);
        fastaName = fastaName.substr(0, fastaName.find_last_of('.'));    
        std::cout << "Parsing " + fastaPath << std::endl;

        OutputBed oBed{outputDir, fastaPath};
        OutputRtr oRtr{outputDir, fastaPath};

        FastaReader fr(fastaPath, 1000);
        for (auto &chrom : *fr.read()) {
            seqFileMap[fastaName + *chrom.first] = outputDir + fastaName + "_" + chrom.first->substr(1);

            // Scoring
            std::cout << "Scoring" << *chrom.first << std::endl;
            ScorerTr st(*chrom.second, 13, 250, 34000);

            // Forward Merging
            std::cout << "Forward Merge" << *chrom.first << std::endl;
            Merger forwardM(st.getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                            LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);

            // Backward Merging
            std::cout << "Backward Merge" << *chrom.first <<  std::endl;
            Merger backwardM(st.getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
                            LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);
            
            //writeStretches(seqFileMap[fastaName + *chrom.first] + "_Forward", *forwardM.getStretchVec());
            //writeStretches(seqFileMap[fastaName + *chrom.first] + "_Backward", *backwardM.getStretchVec());

            // Applying Classifier
            std::cout << "Detecting" << *chrom.first << std::endl;
            Detector dt{red, *chrom.second};
            //writeScores(seqFileMap[fastaName + *chrom.first] + ".scr", red.score(*chrom.second), *chrom.first);
            //writeRepeats(seqFileMap[fastaName + *chrom.first] + ".rpt", dt.getRepeats(), *chrom.first);


            auto fElement = dt.apply(*forwardM.getStretchVec());
            auto bElement = dt.apply(*backwardM.getStretchVec());

            //writeEleFile(outputDir + fastaName + "_" + chrom.first->substr(1) + ".fele", fElement);
            //writeEleFile(outputDir + fastaName + "_" + chrom.first->substr(1) + ".bele", bElement);

            // Matching
            std::cout << "Matching" << *chrom.first << std::endl;
            Matcher mat{fElement, bElement, red, *ic, chrom.second};

            auto rtVecPtr = mat.getRtVec();


            PostProcess post{*rtVecPtr};
            post.apply();

            // Filtering
            // std::cout << "Filtering" << *chrom.first << std::endl;
            // Filter filter(*rtVecPtr, red, chrom.second);
            // filter.apply();


            std::string chromName = chrom.first->substr(1, std::string::npos);
            chromName = chromName.substr(0, chromName.find(" "));
            oBed.write(chromName, *rtVecPtr);
            oRtr.write(chromName, *rtVecPtr);
        }
    }

    delete ic;

}

// Appending elements to database file and writing elements to their own file for later retrieval
int writeToFiles(std::string dbPath, std::string elePath, std::vector<Element> &eleVec, std::string *seq, int fastaId, int id) {
    std::fstream dbFile{dbPath, ios_base::app};
    std::ofstream eleFile{elePath};
    for (auto &ele : eleVec) {
        if (ele.getSize() >= LtrParameters::MIN_LTR) {
            dbFile << ">" << fastaId << "_" << id << std::endl;
            dbFile << seq->substr(ele.getStart(), ele.getEnd() - ele.getStart()) << std::endl;

            eleFile << ele << std::endl;
            id += 1;
        }
    }
    dbFile.close();
    eleFile.close();
    return id;
}

void writeEleFile(std::string elePath, std::vector<Element> &eleVec) {
    std::ofstream eleFile{elePath};
    for (auto &ele : eleVec) {
        eleFile << ele << std::endl;
    }
    eleFile.close();
}

std::vector<Element> readEleFile(std::string elePath) {
    std::vector<Element> r;
    std::fstream eleFile{elePath};
    std::string line;
    std::string extra;
    std::string stretchData;

    while (std::getline(eleFile, line)) {
        std::vector<std::tuple<int, int, int>> stretchDataVec;
        bool isForward = true;
        std::stringstream ss{line};
        ss >> extra;
        while (ss >> stretchData) {
            if (stretchData != "-" && stretchData != "+") {
                int medianHeight;
                int start;
                int end;

                std::stringstream ssStretch{stretchData};
                std::string token;

                std::getline(ssStretch, token, ',');
                medianHeight = std::stoi(token.substr(1));

                std::getline(ssStretch, token, ',');
                start = std::stoi(token);

                std::getline(ssStretch, token);
                end = std::stoi(token.substr(0, token.size() - 1));

                stretchDataVec.push_back(std::tuple<int, int, int>{medianHeight, start, end});
                
            }
            else {
                isForward = stretchData == "+"? true:false;
            }
        }
        std::vector<Stretch*> stretchVec;

        for (auto &x : stretchDataVec) {
            stretchVec.push_back(new Stretch{std::get<1>(x), std::get<2>(x), 0, isForward});
            stretchVec.back()->setMedianHeight(std::get<0>(x));
        }

        Element e{*stretchVec.at(0)};
        int x = 1;
        while (x < stretchVec.size()) {
            e.merge(*stretchVec.at(x));
            x++;
        }

        r.push_back(e);
    }
    return r;

}

void writePrediction(std::string path, std::vector<int> pred) {
    std::ofstream file{path};

    for (auto i : pred) {
        file << i << std::endl;
    }

    file.close();

}

void writeFeatures(std::string path, Matrix features) {
    features.printToFile(path);
}

void writeScores(std::string path, std::vector<int> scoreVec, std::string chrom) {
    int count = 0;
    std::ofstream file{path};
    file << chrom << std::endl;

    for (auto s : scoreVec) {
        if (count == 50) {
            file << std::endl;
            count = 0;
        }
        file << s << " ";
        count++;
    }
    file.close();
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

void writeRepeats(std::string filePath, std::vector<ILocation*>* repeatVec, std::string chrom) {
    std::ofstream file{filePath};
    for (auto r : *repeatVec) {
        file << chrom << "\t" << r->getStart() << "\t" << r->getEnd() + 1 << std::endl;
    }
    file.close();
}
