#include <iostream>
#include <array>
#include <cmath>

/* fun declaration */
/* primitives */
double calc_Cs(double gamma, double p, double rho){
	return sqrt((gamma)*(p/rho));
}

double calc_alpha(double u, double Cs){
	return abs(u) + Cs;
}
/* primitives */
/* array primitives */
bool calc_Cs_CArray(double *gamma, double *p, double *rho, 
		    double* csArray, int arrayLength){
	for (int i=0; i<arrayLength; i++){
		csArray[i] = sqrt((gamma[i]) * (p[i]/rho[i]));
	}
	return true;
}

bool calc_alpha_CArray(double *u, double *Cs, 
		       double *aMaxArray, int arrayLength){
	for (int i=0; i<arrayLength; i++){
		aMaxArray[i] = abs(u[i]) + Cs[i];
	}
	return true;
}

double calc_alphaMax_CArray(double *u, double *Cs,
			    double *aMaxArray, int arrayLength){
	double aMax = 0;
	double aNew, aOld;
	for (int i=0; i<arrayLength-1; i++){
		aOld = calc_alpha(u[i], Cs[i]);
		aNew = calc_alpha(u[i+1], Cs[i+1]);		
		if (aNew > aOld) aMax = aNew;	
	}

	return aMax;
}
/* array primitives */

double calc_deltaT(double CFL, double deltaX, double aMax){
	return (CFL*deltaX)/aMax;
}

/* FORCE Solver
 *
 *	 FORCE	 1    LF    n   n       RI    n   n 
 *	f     =	 - [ f    (q , q   ) + f    (q , q  ) ] 
 *	i+1/2 	 2   i+1/2  i   i+1     i+1/2 i   i+1
 * 
 * calculate Lax-Friedrihcs
 *     
 *	LF     n   n      1     n       n      1  Dx   n    n
 *     f     (q , q  ) =  - [f(q ) + f(q  )] + -  -- (q  - q    )
 *      i+1/2  i   i+1    2     i       i+1    2  Dt   i    i+1
 * */ 
double calc_microLF(double q, double qNext, double fq, double fqNext,
		    double deltaX, double deltaT){	
	double halfDelta = 0.5 * (deltaX/deltaT);
	(0.5 * (fq + fqNext)) +(halfDelta*(q - qNext));
}

double FORCE_LF(std::array<double,3> q, std::array<double,3> qNext,
		std::array<double,3> fq, std::array<double,3>fqNext,
		std::array<double,3> & F_lf){
	F_lf[0] = calc_microLF(q[0], qNext[0], fq[0], fqNext[0], 1, 1);
	F_lf[1] = calc_microLF(q[1], qNext[1], fq[1], fqNext[1], 1, 1);
	F_lf[2] = calc_microLF(q[2], qNext[2], fq[2], fqNext[2], 1, 1);
}

/* calculate Richtmyer */
double calc_microR


/* fun declaration */

int main(void){
	


	return 1;
}


