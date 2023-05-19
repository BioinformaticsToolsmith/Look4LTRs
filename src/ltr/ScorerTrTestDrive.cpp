/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Test.cpp
 *
 *  Created on: Oct 4, 2022
 *      Author: Anthony B. Garza.
 *      Edited/reviewed by Hani Z. Girgis
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

#include "LtrParameters.h"
#include "ScorerTr.h"
#include "../FastaReader.h"


#include <filesystem>
#include <string>

int main(int argc, char *argv[])
{

	std::string input(argv[1]);
	std::string outDir(argv[2]);
	std::string fastaName = filesystem::path(argv[1]).stem().string();

	FastaReader fr(input, 1000);
	std::string *x = fr.read()->at(0).second;

	ScorerTr st(*x, 13, 250, 34000);

	st.printBackwardScores(outDir + "/" + fastaName + ".bscr");
	st.printForwardScores(outDir + "/" + fastaName + ".fscr");
}
