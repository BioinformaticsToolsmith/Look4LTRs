/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Segment.h
 * 
 *  Created on: Sep 23, 2022
 *      Author: Hani Z. Girgis, PhD and Anthony B. Garza.
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
#ifndef SRC_LTR_STRETCH_H_
#define SRC_LTR_STRETCH_H_

#pragma once

#include <iostream>
#include <assert.h>

class Stretch {
	public:

		// Possible values of mark
		static const int K = 0;  // stands for Keep
		static const int D = 1;  // stands for Delete
		static const int I = 2;  // stands for Interrupt

		// Constructor, Destructor, Copy, Move
		Stretch(int start, int end, int mark, bool isForward);
		virtual ~Stretch();
		Stretch(const Stretch &other);
		// Segment(Segment &&other);

		// Given another stretch, calculate the gap between the two stretches
		int calculateGap(const Stretch &anotherStretch) const;

		// Given another stretch and a mark, merge the two stretches and return the new stretch with the given mark
		Stretch merge(Stretch &other, int mark);

		// Builds a stretch that points to to this stretch
		Stretch buildMatch();

		bool isPartOf(int s, int e, double overlap_min) const;

		// Getter and Setters
		int getEnd() const;
		void setEnd(int end);

		int getMark() const;
		void setMark(int mark);

		int getStart() const;
		void setStart(int start);

		int getMedianHeight() const;
		void setMedianHeight(int medianHeight);

		bool getIsForward() const;

		int getSize() const;

	private:
		int start;
		int end;  // end is exclusive
		int mark;
		int medianHeight = 0; // median height of the stretch; distance to matching stretch
		bool isForward;	// is pointing forward

		void check();  // checks for valid arguments
};

// Returns an ostream containing the start, end, and mark of the stretch
std::ostream& operator<<(std::ostream &os, const Stretch &n);

#endif /* SRC_LTR_SEGMENT_H_ */
