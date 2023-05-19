#include "../ltr/LtrUtility.h"
#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>

std::string reverseComplement(std::string &seq) {
    std::unordered_map<char, char> complementTable = {
        {'A', 'T'},
        {'T', 'A'},
        {'C', 'G'},
        {'G', 'C'},
        {'N', 'N'}
    };
    std::string r = seq;
    std::reverse(r.begin(), r.end());
    std::transform(r.begin(), r.end(), r.begin(), [&complementTable](char c) {
        return complementTable.at(c);
    });
    return r;
}

int main() {

    std::string a = "ATCG";
    std::cout << reverseComplement(a) << std::endl;

}