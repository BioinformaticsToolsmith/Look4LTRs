#include "StretchFeature.h"
#include "Stretch.h"
#include "../Matrix.h"

#include <vector>

Stretch makeStretch(int start, int end, int height, int mark = 0);

int main() {
    std::vector<Stretch> stretchVec = {
        
         makeStretch(2, 11, 700),

         makeStretch(13, 26, 60),
         makeStretch(26, 35, 66),

         makeStretch(40, 54, 900),

         makeStretch(54, 110, 60), 


         makeStretch(200, 275, 125),

         makeStretch(325, 350, 125),
         makeStretch(355, 400, 130),

         
         makeStretch(450, 475, 300),
         makeStretch(500, 525, 325)
         };

    std::string bedPath = "/home/anthony/Documents/LtrDetector_v2/c++/Oct6/TestCases/Bed/test.bed";

    StretchFeature a(LtrParameters::MEAN_VECTOR, LtrParameters::STD_VECTOR);
    Matrix features = a.buildMatrix(stretchVec);
    features = a.scale(features);

    return 0;
}

Stretch makeStretch(int start, int end, int height, int mark) {
    Stretch a = {start, end, mark};
    a.setMedianHeight(height);
    return a;
}
