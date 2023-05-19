/*
 * Test2.cpp
 *
 *  Created on: Nov 25, 2020
 *      Author: Hani Z. Girgis, PhD
 */
#include <iostream>
#include <math.h>

double gKernel(double id) {
	double dist = 1.0 - id;
	return exp(-dist * dist);
}

int main(int argc, char *argv[]) {
	for (double i = 0; i < 1.0; i += 0.1) {
		std::cout << i << " " << gKernel(i) << std::endl;
	}
}

