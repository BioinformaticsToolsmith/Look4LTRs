#include <vector>
#include <algorithm>
#include <iostream>

int median(std::vector<int> &scoreVec) {
    int r;
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
int main() {

    // 1 2 3 4 5 6 7 8 9 10
    std::vector<int> a{5, 1, 0, 0, 0, 0, 8, 9, 0, 0, 2, 4, 7, 0, 0, 3, 6, 10};

    std::vector<int> b{5, 0, 0, 0, 0, 1, 8, 9, 0, 0, 2, 4, 7, 3, 0, 0, 6};

    std::cout << median(a) << std::endl << median(b) << std::endl;

}

