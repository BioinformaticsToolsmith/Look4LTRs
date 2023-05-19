
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <memory>

int main() {

    std::ifstream j;
    j.open("/home/transposons/Projects/Identity/cases/redmulti/test.scr");
    std::string fasta = "test";
    std::unordered_map<std::string, std::shared_ptr<std::vector<int>> > m;

    if (j.is_open()) {

        std::vector<int> scoreVec;
        std::string seqName;
        std::string line;
        bool isSeqStarted = false;
        while (std::getline(j, line)){
            std::cout << line << std::endl;
            if (line.size() == 0 && isSeqStarted == true) {
                m[fasta + ":" + seqName] = std::make_shared<std::vector<int>>(std::vector<int>{scoreVec.begin(), scoreVec.end()});
                
                // Don't reallocate memory; may need to fill up again
                scoreVec.clear();

                isSeqStarted = false;
            }
            else if (line.at(0) == '>') {
                seqName = line;
                isSeqStarted = true;
            } 
            else {
                std::stringstream s(line);
                int score;
                while (s >> score) {
                    scoreVec.push_back(score);
                }
            }
        }
    }
    else {
        throw std::exception();
    }

}