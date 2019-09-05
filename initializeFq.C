#include "header.h"
#include <cmath>
#include <iostream>

//double divide(double, double);

void initializeFq(std::vector<std::array<double, 3>> &fq, int T){
/*
 * Function for f(q) initialization. 
 * It works with the arrays for r, u and p that have been passed as constants. 
 * In order to work with new data it should be modified to accept the three arrays as arguments. 
 * */
	
	//const double constants::gamma;
	//const std::array<double,10> r = {1, 1, 1, 1, 5.99924, 0.125, 1, 1, 1, 5.99242};
	//const std::array<double,10> u = {0, -2, 0 ,0, 19.5975, 0, 2, 0, 0, -6.19633};
	//const std::array<double,10> p = {1, 0.4, 1000, 0.01, 460.894, 0.1, 0.4, 0.01, 100, 46.095};
 
	double E, t, pt, rt, ut;       // E is energy, t is for time.
        int n = fq.size();   // number of cells.
	
        for (int i = 0; i < n; i++){
        	if (i < n/2){
                	t = T-1;
                }
                else{
                        t = T+4;
                }
		pt = constants::p[t];
		rt = constants::r[t];
		ut = constants::u[t];
                E = pt/constants::gamma + ((1/2)*rt*pow(ut, 2)); 
		// Compute E

                fq.at(i)[0] = rt * ut;
                fq.at(i)[1] = (rt * pow(ut,2)) + pt;
		fq.at(i)[2] = (E + pt)*ut;
//		std::cout << fq.at(i)[0] << "\t" << fq.at(i)[1] << "\t" << fq.at(i)[2] << "\n";

       }
}

/*double divide(double a, double b){
	return double(a/b);
}*/
