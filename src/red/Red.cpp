/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Red
 *
 *  Created on: Dec 26, 2022
 *      Author: Hani Z. Girgis, Anthony B. Garza
 * Reviewer: @@@@@@@@@@@@@@@@@@@@@@@@ Need to be reviewed @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *   Purpose:
 *
 *
 * Academic use: Affero General Public License version 1.
 *
 * Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 *
 * Copyright (C) 2022 by the authors.
 */

#pragma once

#include "Red.h"

Red::Red(std::string _gnm, std::optional<int> _cor, std::optional<int> _k,
         std::optional<int> _ord, std::optional<double> _gau, std::optional<double> _thr,
         std::optional<int> _min, std::optional<std::string> _cnd)
{
    assert(std::filesystem::exists(_gnm));
    assert(_cor.has_value() ? _cor.value() >= 1 : true);
    assert(_thr.has_value() ? _thr.value() >= 1 : true);
    assert(_min.has_value() ? _min.value() >= 0 : true);
    if (_cnd.has_value())
    {
        assert(std::filesystem::exists(_cnd.value()));
    }

    gnm = _gnm;
    cor = _cor.has_value() ? _cor.value() : Util::CORE_NUM;
    thr = _thr.has_value() ? _thr.value() : 2.0;
    thr = static_cast<int>(thr) == 1 ? 1.5 : thr;
    min = _min.has_value() ? _min.value() : 3;

    // Deal with the value of k
    if (_k.has_value())
    {
        k = _k.value();
    }
    else
    {
        calcGnmLen();
    }
    if (k > 15)
    {
        std::cout << "Due to a memory constraint, k is set to 15.";
        std::cout << std::endl;
        k = 15;
    }
    if (k < 12)
    {
        std::cout << "Due to a statistical consideration, k is set to 12.";
        std::cout << std::endl;
        k = 12;
    }
    std::cout << std::endl;

    // Deal with the value of order
    if (_ord.has_value())
    {
        ord = _ord.value();
        assert(ord >= 1);
    }
    else
    {
        calcMarkovOrd();
    }

    // Deal with Gaussian Window
    if (_gau.has_value())
    {
        gau = _gau.value();
        if (gau < 20)
        {
            gau = 20;
            std::cout << "Enlarging the smoothing window to be large enough." << std::endl;
        }
        else if (gau > 40)
        {
            gau = 40;
            std::cout << "Decreasing the smoothing window to be small enough." << std::endl;
        }
    }
    else
    {
        calcGauWidth();
    }

    std::cout << "Red will run with the following parameters:" << std::endl
              << "\tgnm: " << gnm << std::endl
              << "\tcor: " << cor << std::endl
              << "\tlen: " << k << std::endl
              << "\tord: " << ord << std::endl
              << "\tgau: " << gau << std::endl
              << "\tthr: " << thr << std::endl
              << "\tmin: " << min << std::endl;
    if (_cnd.has_value())
    {
        std::cout << "\tcnd: " << _cnd.value() << std::endl;
    }

    // Train
    if (_cnd.has_value())
    {
        trainer = new Trainer(gnm, ord, k, gau, thr, _cnd.value(), min);
    }
    else
    {
        trainer = new Trainer(gnm, ord, k, gau, thr, min);
    }
}

Red::~Red()
{
    delete trainer;
}

void Red::calcGnmLen()
{
    // Calculate the size of the genome
    long genomeLength = 0;

    std::vector<std::string> *fileList = new std::vector<std::string>();
    Util::readChromList(gnm, fileList, "fa");
    std::cout << "Calculating the length, k, of the k-mer ";
    std::cout << "based on the input genome ... " << std::endl;

    #pragma omp parallel for schedule(dynamic) num_threads(Util::CORE_NUM) reduction(+: genomeLength)
    for (unsigned int i = 0; i < fileList->size(); i++)
    {

        #pragma omp critical
        {
            std::cout << "\t" << fileList->at(i) << std::endl;
        }
        ChromListMaker *maker = new ChromListMaker(fileList->at(i));
        const std::vector<Chromosome *> *chromList = maker->makeChromList();
        for (unsigned int h = 0; h < chromList->size(); h++)
        {
            genomeLength += chromList->at(h)->getEffectiveSize();
        }
        delete maker;
    }
    fileList->clear();
    delete fileList;

    double temp = log(genomeLength) / log(4.0);

    k = floor(temp);
    std::cout << "The recommended k is " << k << "." << std::endl;
}

void Red::calcMarkovOrd()
{
    ord = floor(k / 2.0) - 1;

    std::cout << "Using the default background order: " << ord << ".";
    std::cout << std::endl;
}

void Red::calcGauWidth()
{
    std::cout << "Calculating GC content ..." << std::endl;

    // 1: Count the gc content of the input genome
    long genomeLength = 0;
    long genomeGc = 0;
    std::vector<std::string> *fileList = new std::vector<std::string>();
    Util::readChromList(gnm, fileList, "fa");

#pragma omp parallel for num_threads(Util::CORE_NUM) schedule(dynamic) reduction(+ \
                                                                                 : genomeGc, genomeLength)
    for (unsigned int i = 0; i < fileList->size(); i++)
    {
#pragma omp critical
        {
            std::cout << "\t" << fileList->at(i) << std::endl;
        }

        ChromListMaker *maker = new ChromListMaker(fileList->at(i));
        const std::vector<Chromosome *> *chromList = maker->makeChromList();
        for (unsigned int h = 0; h < chromList->size(); h++)
        {
            genomeGc += chromList->at(h)->getGcContent();
            genomeLength += chromList->at(h)->getEffectiveSize();
        }

        delete maker;
    }
    fileList->clear();
    delete fileList;

    // 2: Calculate the gc content of the input genome
    double gc = 100.00 * genomeGc / genomeLength;

    gau = gc < 33 || gc > 67 ? 40.0 : 20.0;

    std::cout << "Using the default half width: " << gau;
    std::cout << " based on the GC content of " << gc << std::endl;
}

void Red::printTable(std::string fileName)
{
    assert(std::filesystem::exists(fileName));
    trainer->printTable(fileName);
}

void Red::printModel(std::string fileName)
{
    assert(std::filesystem::exists(fileName));
    trainer->printHmm(fileName);
}

// @@@@@@@@@@@ TODO: MAKE IT PARALLEL @@@@@@@@@@@@@@@@@@@@@@@@@
const std::vector<ILocation*>* Red::predictRepeats(std::string &seq)
{
    ChromosomeOneDigit chrom(seq, "temp");
    
    HMM *copyHMM = new HMM(*trainer->getHmm());

    // Scan the forward strand
    Scanner *scanner = new Scanner(copyHMM, k, &chrom,
                                    trainer->getTable());

    // Scan the reverse complement
    chrom.makeRC();
    Scanner *scannerRC = new Scanner(copyHMM, k, &chrom,
                                        trainer->getTable());
    scannerRC->makeForwardCoordinates();
    scanner->mergeWithOtherRegions(scannerRC->getRegionList());
    delete scannerRC;
    chrom.makeRC();

    // Scan the reverse
    chrom.makeR();
    Scanner *scannerR = new Scanner(copyHMM, k, &chrom,
                                    trainer->getTable());
    scannerR->makeForwardCoordinates();
    scanner->mergeWithOtherRegions(scannerR->getRegionList());
    delete scannerR;

    auto regionList = scanner->getRegionList();
    std::vector<ILocation*>* r = new std::vector<ILocation*>{regionList->size(), nullptr};
    for (int i = 0; i < regionList->size(); i++) {
        r->at(i) = new Location(* regionList->at(i) );
    }

    delete scanner;  
    delete copyHMM; 

    return r;  

}

void Red::scan(std::string repeatPath, int format, std::optional<std::string> scorePath,
          std::optional<std::string> maskPath, std::optional<std::string> dirPath)
{

    // pre-conditons
    assert(std::filesystem::exists(repeatPath));
    assert(format == 1 || format == 2);

    if (scorePath.has_value())
    {
        assert(std::filesystem::exists(scorePath.value()));
    }
    if (maskPath.has_value())
    {
        assert(std::filesystem::exists(maskPath.value()));
    }
    if (dirPath.has_value())
    {
        assert(std::filesystem::exists(dirPath.value()));
    }

    // Stage 4: Scan
    std::cout << std::endl
              << std::endl;
    std::cout << "Stage 4: Scanning ..." << std::endl;
    std::vector<std::string> *fileList = new std::vector<std::string>();
    Util::readChromList(gnm, fileList, std::string("fa"));
    if (dirPath.has_value())
    {
        Util::readChromList(dirPath.value(), fileList, std::string("fa"));
    }

    unsigned long genomeLen = 0;
    unsigned long repeatLen = 0;

    unsigned int chromCount = fileList->size();
#pragma omp parallel for schedule(dynamic) num_threads(cor)
    for (unsigned int i = 0; i < chromCount; i++)
    {
#pragma omp critical
        {
            cout << "Scanning: " << fileList->at(i) << endl;
        }
        // Output file name
        std::string path(fileList->at(i));
        int slashLastIndex = path.find_last_of(Util::fileSeparator);
        int dotLastIndex = path.find_last_of(".");
        std::string nickName = path.substr(slashLastIndex + 1,
                                      dotLastIndex - slashLastIndex - 1);

        // Process each sequence with the ith file
        ChromListMaker *maker = new ChromListMaker(fileList->at(i));
        const std::vector<Chromosome *> *chromList =
            maker->makeChromOneDigitList();

        ChromListMaker *oMaker = new ChromListMaker(fileList->at(i));
        const std::vector<Chromosome *> *oChromList;
        if (maskPath.has_value())
        {
            oChromList = oMaker->makeChromList();
        }

        for (unsigned int h = 0; h < chromList->size(); h++)
        {
            ChromosomeOneDigit *chrom =
                dynamic_cast<ChromosomeOneDigit *>(chromList->at(h));
            genomeLen += chrom->size();

            HMM *copyHMM = new HMM(*trainer->getHmm());

            // Scan the forward strand
            Scanner *scanner = new Scanner(copyHMM, k, chrom,
                                           trainer->getTable());

            // Scan the reverse complement
            chrom->makeRC();
            Scanner *scannerRC = new Scanner(copyHMM, k, chrom,
                                             trainer->getTable());
            scannerRC->makeForwardCoordinates();
            scanner->mergeWithOtherRegions(scannerRC->getRegionList());
            delete scannerRC;
            chrom->makeRC();

            // Scan the reverse
            chrom->makeR();
            Scanner *scannerR = new Scanner(copyHMM, k, chrom,
                                            trainer->getTable());
            scannerR->makeForwardCoordinates();
            scanner->mergeWithOtherRegions(scannerR->getRegionList());
            delete scannerR;

            repeatLen += scanner->getTotalRegionLength();

            //@@ The chromosome now has the sequence of the reverse strand
            // The actual strand is calculated if the user requested the scores.

            // Print according to the user's requests
            bool canAppend = (h == 0) ? false : true;

            if (scorePath.has_value())
            {
                // Calculate the forward strand from the reverse
                chrom->makeR();

                std::string scoFile = scorePath.value() + Util::fileSeparator + nickName + ".scr";
                if (!canAppend)
                {
                    #pragma omp critical
                    {
                        cout << "Printing scores to: " << scoFile << endl;
                    }
                }
                // Make sure to print the original E-values not their logarithm
                Scorer *scorer = new Scorer(chrom, trainer->getTable());
                scorer->printScores(scoFile, canAppend);
                delete scorer;
            }

            std::string ext(".rpt");
            if (format == 2)
            {
                ext = std::string(".bed");
            }
            std::string rptFile = repeatPath + Util::fileSeparator + nickName + ext;
            if (!canAppend)
            {
                #pragma omp critical
                {
                    cout << "Printing locations to: " << rptFile << endl;
                }
            }
            scanner->printIndex(rptFile, canAppend, format);

            if (maskPath.has_value())
            {
                std::string mskFile = maskPath.value() + Util::fileSeparator + nickName + ".msk";
                if (!canAppend)
                {
                    #pragma omp critical
                    {
                        cout << "Printing masked sequence to: " << mskFile << endl;
                    }
                }
                Chromosome *oChrom = oChromList->at(h);
                scanner->printMasked(mskFile, *oChrom, canAppend);
            }

            // Free memory
            delete scanner;
            delete copyHMM;
        }

        delete maker;
        delete oMaker;
    }

    cout << "Genome length: " << genomeLen;
    cout << " - Repeats length: " << repeatLen;
    cout << " - Repeats content: "
         << 100.00 * ((double)repeatLen / genomeLen) << "%" << endl;

    // Free memory
    fileList->clear();
    delete fileList;
}

std::vector<int> Red::score(std::string &seq){
    ChromosomeOneDigit chrom(seq, "temp");
    Scorer scorer(&chrom, trainer->getTable());
    auto scoreVec = scorer.getScores();
    std::vector<int> r(scoreVec->begin(), scoreVec->end());
    return r;
}

double Red::calcPercent(std::string &seq) {
    auto scoreVec = score(seq);
    double count = 0.0;
    for (auto x : scoreVec) {
        if (x != 0) {
            count++;
        }
    }

    return count / scoreVec.size();

}

int Red::calcMedianScore(std::string &seq) {
    int r;
    auto scoreVec = score(seq);
    std::vector<int> noZeroVec;
    for (auto x : scoreVec) {
        if (x != 0) {
            noZeroVec.push_back(x);
        }
    }

    int middle = noZeroVec.size() / 2;
    if (noZeroVec.size() % 2 == 1) {
        std::nth_element(noZeroVec.begin(), noZeroVec.begin() + middle, noZeroVec.end());
        r = noZeroVec.at(noZeroVec.size() / 2);
    }
    else {
        std::nth_element(noZeroVec.begin(), noZeroVec.begin() + middle, noZeroVec.end());
        std::nth_element(noZeroVec.begin(), noZeroVec.begin() + middle - 1, noZeroVec.end());
        r = static_cast<int> ((noZeroVec.at(middle) + noZeroVec.at(middle - 1)) / 2.0);
    }
    return r;
}

double Red::calcMeanScore(std::string &seq, bool withZero) {
    double r;
    auto scoreVec = score(seq);
    auto *scorePtr = &scoreVec;
    int sum = 0;
    std::vector<int> noZeroVec;
    if (withZero == false) {
        for (auto x : scoreVec) {
            if (x != 0) {
                noZeroVec.push_back(x);
            }
        }       
        scorePtr = &noZeroVec;
    }
    
    for (auto x : *scorePtr) {
        sum += x;
    }
    r = sum / double(scorePtr->size());
    return r;
}


int Red::getK(){
    return k;
}
