#include "LtrParameters.h"
#include "Merger.h"
#include "ScorerTr.h"
#include "../FastaReader.h"

#include <filesystem>
#include <string>

std::vector<int> fillScores(std::vector<std::tuple<int, int, int>> stretches, int remainingGap = 0);
std::vector<int> readScores(std::string fileName);

int main(int argc, char *argv[])
{

	std::string input(argv[1]);
	std::string outDir(argv[2]);
	std::string fastaName = filesystem::path(argv[1]).stem().string();

	FastaReader fr(input, 1000);
	std::string *x = fr.read()->at(0).second;

	ScorerTr st(*x, 13, 250, 34000);

    Merger mf(st.getForwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
        LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, true);
    Merger mb(st.getBackwardScores(), LtrParameters::MIN_STRETCH, LtrParameters::MAX_GAP,
        LtrParameters::SIM_MARGIN, LtrParameters::INT_MARGIN, false);

    mf.printScores(outDir + "/" + fastaName + ".fscr");
    mb.printScores(outDir + "/" + fastaName + ".bscr");

return 0;
}

    // // Contains four stretches, all mergable with each other
	// std::vector<std::tuple<int, int, int>> mergable{
    //     std::make_tuple(30, 90, 2000),
    //     std::make_tuple(100, 110, 2003),
    //     std::make_tuple(115, 120, 1999),
    //     std::make_tuple(125, 150, 2000),
    // };
	// // The correct answer is 30 150 2000

    // // Contains five stretches, all not mergable with each other
	// std::vector<std::tuple<int, int, int>> unmergable{
    //     std::make_tuple(50, 80, 2000),
    //     std::make_tuple(90, 100, 3000),
    //     std::make_tuple(130, 140, 1000),
    //     std::make_tuple(180, 200, 750),
    //     std::make_tuple(220, 280, 1500)
    // };

    // /*
    //  * Contains three stretches with the middle being an interrupt
    //  * and the others being keep segments
    //  */
    // std::vector<std::tuple<int, int, int>> interruptMerge{
    //     std::make_tuple(30, 90, 2000),
    //     std::make_tuple(95, 100, 400),
    //     std::make_tuple(100, 150, 2010)
    // };

    // /*
    //  * Contains three stretches with the middle being an interrupt,
    //  * the first being a keep, and the last being a delete
    //  */
    // std::vector<std::tuple<int, int, int>> interruptMergeF{
    //     std::make_tuple(30, 90, 2000),
    //     std::make_tuple(95, 100, 400),
    //     std::make_tuple(100, 105, 2010)
    // };

    // /*
    //  * Contains three stretches with the middle being an interrupt,
    //  * the first being a keep, and the last being a delete
    //  */
    // std::vector<std::tuple<int, int, int>> interruptMergeB{
    //     std::make_tuple(85, 90, 2000),
    //     std::make_tuple(95, 100, 400),
    //     std::make_tuple(100, 150, 2010)
    // };

    // /*
    //  * Contains three stretches with similar scores separated
    //  * by large gaps that prevent mergin
    //  */
    // std::vector<std::tuple<int, int, int>> gapUnmergable{
    //     std::make_tuple(30, 100, 2000),
    //     std::make_tuple(200, 250, 2002),
    //     std::make_tuple(340, 380, 2004)
    // };

    // /*
    //  * Contains two similar keep stretches separated by two
    //  * non-similar delete stretches
    //  */
    // std::vector<std::tuple<int, int, int>> twoInterruptMergable{
    //     std::make_tuple(30, 100, 2000),
    //     std::make_tuple(110, 115, 3000),
    //     std::make_tuple(120, 125, 1000),
    //     std::make_tuple(130, 180, 2001)
    // };

    // /*
    //  * Contains two similar stretches separated by two
    //  * non-similar delete stretches that can only be merged
    //  * if the two inner delete stretches are marked as interrupts
    //  * and are deleted. Furthermore, the merge must be forward
    //  */
    // std::vector<std::tuple<int, int, int>> twoInterruptUnmergableF{
    //     std::make_tuple(30, 100, 2000),
    //     std::make_tuple(110, 115, 3000),
    //     std::make_tuple(120, 125, 1000),
    //     std::make_tuple(130, 155, 2001)
    // };


    // /*
    //  * Contains two similar stretches separated by two
    //  * non-similar delete stretches that can only be merged
    //  * if the two inner delete stretches are marked as interrupts
    //  * and are deleted. Furthermore, the merge must be backwards
    //  */
    // std::vector<std::tuple<int, int, int>> twoInterruptUnmergableB{
    //     std::make_tuple(90, 100, 2000),
    //     std::make_tuple(110, 115, 3000),
    //     std::make_tuple(120, 125, 1000),
    //     std::make_tuple(130, 190, 2001)
    // };

    // std::vector<int> stretches = fillScores(interruptMerge, 10);
    // std::vector<int> stretches = readScores("/home/transposons/Projects/Identity/src/ltr/scores_forward.txt");


std::vector<int> readScores(std::string fileName) {
	std::ifstream in(fileName.c_str());
	if (!in.good()) {
		std::cerr << "Could not open " << fileName << std::endl;
		throw std::exception();
	}
	int score;
	std::vector<int> v;
	while (in >> score) {
		//std::cout << score << std::endl;
		v.push_back(score);
	}

	in.close();
	return v;
}

std::vector<int> fillScores(std::vector<std::tuple<int, int, int>> stretches, int remainingGap) {

	std::vector<int> scoreVec(std::get<1>(stretches.back()) + remainingGap, 0);
	int score = 0;
	for (int x = 0; x < stretches.size(); x++) {
		score = std::get<2>(stretches[x]);
		for (int i = std::get<0>(stretches[x]); i < std::get<1>(stretches[x]);
				i++) {
			scoreVec[i] = score;
		}
	}

	// for (int x = 0; x < scoreVec.size(); x++) {
	// 	std::cout << x << " " << scoreVec[x] << std::endl;
	// }

	return scoreVec;
}
