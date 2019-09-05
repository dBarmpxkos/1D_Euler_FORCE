#include "header.h"
#include <cmath>
#include <iostream>


void initializeq(std::vector<std::array<double, 3>> &q, int T){
/*
 * Function for q initialization. 
 * It works with the arrays for r, u and p that have been passed as constants. 
 * In order to work with new data it should be modified to accept the three arrays as arguments. 
 * */
   
		double E, t;       // E is energy.
                int n = q.size();   // number of cells.
                for (int i = 0; i < n; i++){
                        if (i < n/2){
                                t = T-1;
                        }
                        else{
                                t = T+4;
                        }

                        E = (constants::p[t]/constants::gamma) + ((1/2)*constants::r[t]*pow(constants::u[t], 2)); // Compute E

                        q.at(i)[0] = constants::r[t];
                        q.at(i)[1] = constants::r[t]*constants::u[t];
                        q.at(i)[2] = E;
			std::cout << q.at(i)[0] << "\t" << q.at(i)[1] << "\t" << q.at(i)[2] << "\n";
                }
}
