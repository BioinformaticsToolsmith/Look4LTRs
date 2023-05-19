/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * Stretch.cpp
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
 *
 *
 * ToDo: Fix the move constructor
 *
 *
 */
#include "Stretch.h"

Stretch::Stretch(int start, int end, int mark, bool isForward) {

	this->start = start;
	this->end = end;
	this->mark = mark;
	this->isForward = isForward;


	if (!(mark == K || mark == D || mark == I)) {
		std::string msg = "Cons 1: " + std::to_string(mark);
		throw std::invalid_argument(msg);
	}

	check();

}

Stretch::~Stretch() {
	// TODO Auto-generated destructor stub
}

Stretch::Stretch(const Stretch &other) {
	this->start = other.getStart();
	this->end = other.getEnd();
	this->mark = other.getMark();
	this->isForward = other.getIsForward();
	this->medianHeight = other.getMedianHeight();

	if (!(mark == K || mark == D || mark == I)) {
		std::string msg = "Cons 2: " + std::to_string(mark);
		throw std::invalid_argument(msg);
	}

	check();
}

/**
 * Given another stretch, calculate the gap between the two stretches
 * returns: the gap as an integer
 */
int Stretch::calculateGap(const Stretch &other) const {
	int result = 0;
	if (this->start < other.getStart()) {
		result = other.getStart() - this->end;
	}
	else {
		result = this->start - other.getEnd();
	}

	return result;
}

/**
 * Given another stretch and a mark, merge the two stretches
 * returns: a Stretch with the merged coordinates and given mark
 */
Stretch Stretch::merge(Stretch &other, int mark) {
	// if (this->start < other.getStart()) {
	// 	return Stretch{this->start, other.getEnd(), mark, this->isForward};
	// }
	// else {
	// 	return Stretch{other.getStart(), this->end, mark, this->isForward};
	// }

	if (!(mark == K || mark == D || mark == I)) {
		std::string msg = "merge: " + std::to_string(mark);
		throw std::invalid_argument(msg);
	}
	
	return Stretch{std::min(this->start, other.getStart()), std::max(this->end, other.getEnd()), mark, this->isForward};
}

Stretch Stretch::buildMatch() {
	int matchStart = start;
	int matchEnd = end;
	if (isForward) {
		matchStart += medianHeight;
		matchEnd += medianHeight;
	}
	else {
		matchStart -= medianHeight;
		matchEnd -= medianHeight;
	}

	Stretch r(matchStart, matchEnd, mark, !isForward);
	r.setMedianHeight(medianHeight);
	return r;
}


bool Stretch::isPartOf(int s, int e, double overlapMin = 0.50) const {
	assert (s < e);

	bool con_1 = start > s && start < e;
	bool con_2 = end > s && end < e;
	double perc = 0.0;

	if (con_1 || con_2) {
		double overlap = std::min(end, e) - std::max(start, s);
		perc = overlap / double((end - start));
	}

	assert(perc >= 0.0 && perc <= 1.0);
	return perc >= overlapMin;
}

// Stretch::Stretch(Stretch &&other) {
// 	this->start = other.getStart();
// 	this->end = other.getEnd();
// 	this->mark = other.getMark();

// 	check();
// }

/**
 * Check for valid attributes in the object;
 * throws exception if invalid.
 */
void Stretch::check() {
	assert(start < end);

	if (!(mark == K || mark == D || mark == I)) {
		std::string msg = "Invalid mark: " + std::to_string(mark);
		throw std::invalid_argument(msg);
	}

	assert((mark == K || mark == D || mark == I));
}

int Stretch::getEnd() const {
	return end;
}

void Stretch::setEnd(int end) {
	this->end = end;
}

int Stretch::getMark() const {
	return mark;
}

void Stretch::setMark(int mark) {
	if (!(mark == K || mark == D || mark == I)) {
		std::string msg = "setMark: " + std::to_string(mark);
		throw std::invalid_argument(msg);
	}

	this->mark = mark;
	check();
}

int Stretch::getStart() const {
	return start;
}

void Stretch::setStart(int start) {
	this->start = start;
}

std::ostream& operator<<(std::ostream &os, const Stretch &n) {
	os << n.getStart() << " " << n.getEnd() << " " << n.getMark();

	return os;
}


int Stretch::getMedianHeight() const {
	return medianHeight;
}

void Stretch::setMedianHeight(int medianHeight) {
	assert(medianHeight > 0);
	this->medianHeight = medianHeight;
}

bool Stretch::getIsForward() const {
	return isForward;
}


int Stretch::getSize() const {
	return end - start;
}

