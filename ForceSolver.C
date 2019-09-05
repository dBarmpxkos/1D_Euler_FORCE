#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include "header.h"

//void initializeq(std::vector<std::array<double, 3>>, int);
//void initializefq(std::vector<std::array<double, 3>>, int);

//double f (double y);

int main(void){
	//std::cout << constants::gamma;
	int T = -1, cellNo = 1; // Test number and number of cells
	
	while (T < 1 || T > 5){
        	std::cout << "Enter test number (must be a number between 1 and 5): ";
        	std::cin >> T;
	}
	
	while (cellNo % 2 == 1) {
        	std::cout << "Enter number of cells (please, give an even number): ";
        	std::cin >> cellNo; //must be an even number
	}
	

	std::vector<std::array<double, 3>> q(cellNo+2);
	std::vector<std::array<double, 3>> fq(cellNo+1);

	initializeq(q, T);
	initializeFq(fq, T);

}

