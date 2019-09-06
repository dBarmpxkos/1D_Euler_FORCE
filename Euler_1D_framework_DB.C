/* This program is designed to simulate the Euler equations
   in 1 dimension, using the FORCE solver.

   It provides output in ASCII format, suitable for plotting in gnuplot.  
 */

#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <limits>

// Convert from conserved variables to primitive variables
std::array<double, 3> primitive(std::array<double, 3>& q_i)
{
  std::array<double, 3> w_i;

  // TODO
  // Convert the conserved state variables q_i to the primitive variables w_i
  
  return w_i;
}

// Convert primitive state vector to conservative state vector
std::array<double, 3> conservative(std::array<double, 3>& w_i)
{
  std::array<double, 3> q_i;

  // TODO
  // Convert the primitive variables w_i to the conserved variables q_i

  return q_i;
}

// Fill the conservative state vector with the initial data
void initialiseData(std::vector<std::array<double, 3> >& q, 
					int test, double finalT)
{
  int n = q.size();

  if(test == 1)
  {
    for(int i = 0; i < n; ++i)
    {
      std::array<double, 3> w;
      if(i < n/2)
      {
		w = {1.0,0.0,1.0};
      }
      else {
		w = {0.125,0.0,0.1};
      }
      
      std::array<double, 3> q_i = conservative(w);
      q[i][0] = q_i[0];
      q[i][1] = q_i[1];
      q[i][2] = q_i[2];
    }
  }
  else if(test == 2)
  {
    // TODO
    // Insert initial data for this test
  }
  else if(test == 3)
  {
    // TODO
    // Insert initial data for this test
  }
  else if(test == 4)
  {
    // TODO
    // Insert initial data for this test
  }
  else if(test == 5)
  {
    // TODO
    // Insert initial data for this test
  }
  else
  {
    // TODO
    // Add suitable error handling
  }
}
// Compute flux-vector corresponding to given state vector
std::array<double, 3> flux(std::vector<std::array<double, 3>> &fq, int T)
{
  std::array<double, 3> f;
	/* E: energy, t: time */
	double E, t, pt, rt, ut;
	/* Cell number */
    int n = fq.size();   

    for (int i = 0; i < n; i++) {
    	if (i < n/2) {
    		t = T-1;
    	}
    	else {
    		t = T+4;
    	}

    	pt = constants::p[t];
    	rt = constants::r[t];
    	ut = constants::u[t];
    	/* Compute E */
    	E = pt/constants::gamma + ((1/2)*rt*pow(ut, 2)); 

    	fq.at(i)[0] = rt * ut;
    	fq.at(i)[1] = (rt * pow(ut,2)) + pt;
    	fq.at(i)[2] = (E + pt)*ut;
//		std::cout << fq.at(i)[0] << "\t" << fq.at(i)[1] << "\t" << fq.at(i)[2] << "\n";
 
  return f;
}

// Compute the maximum stable time-step for the given data
double computeTimestep(std::vector<std::array<double, 3> >& q, double dx){

  double dt;
  double alphaOld, alphaMax;
  int n = q.size();
  std::array<double, n> Cs;

  // Compute the maximum wavespeed over the entire domain, and use this to compute the timestep
  for(int i=0 ; i < n ; i++)
  {
  	/* q[i][0] and q[i][1] are wrong, i need to find
  	 *  expression for p/rho from q values */
  	Cs[i] = sqrt( gamma * (q[i][0]/q[i][1]) );
  	alphaOld = abs(q[2][i]) + Cs[i];
  	if (alphaOld > alphaMax) 
  		alphaOld = alphaMax;
  }
  dt = (CFL*dx)/alphaMax;

  return dt;
}

// Compute the FORCE flux between two states uL and uR in coordinate direction coord.
std::array<double, 3> FORCEflux(std::array<double, 3>& q_iMinus1, 
								std::array<double, 3>& q_i, 
								double dx, double dt){

	double halfDelta = 0.5 * (dx/dt);
	std::array<double, 3> fluxRM;
	/* Compute the Richtmyer flux using q_i and q_iMinus1 
	   how the f*** only with q's, don't we need fluxes as well? */
	fluxRM = ???

	std::array<double, 3> fluxLF;
	/* Compute the Lax-Friedrichs flux using q_i and q_iMinus1
	   how the f*** only with q's, don't we need fluxes as well? */

	std::array<double, 3> fluxForce;
	// TODO
	// Compute the FORCE flux

	return fluxForce;
}

// Compute the array of fluxes from the given data array
void computeFluxes(std::vector<std::array<double, 3> >& q, std::vector<std::array<double, 3> >& flux, double dx, double dt)
{
  int n = q.size();
  // TODO - consider why this is the choice of index range for the loop
  for(unsigned int i=1 ; i < n ; i++)
  {
    flux[i] = FORCEflux(q[i-1], q[i], dx, dt);
  }
}


int main(void)
{
  // User-defined parameters
  int cells;
  int test;
  // TODO read in cells, and also the desired test
  
  double CFL = 0.9;
  double finalT;
  
  // Set initial vectors:
  std::vector<std::array<double, 3> > q(cells+2);
  std::vector<std::array<double, 3> > flux(cells+2);
  // TODO - consider the extent of the vectors, the flux vector actually contains one more cell than it needs to
  
  initialiseData(q, test, finalT);
  
  // Do simulation
  double t = 0;
  double dx = 1.0 / cells;
  // TODO
  // You may wish to change the 1.0 to a specification of your own domain size
  
  while( t < finalT )
  {
    // TODO
    // Implement this function, and then uncomment
    // computeDomainBoundaries(q)
    
    double dt = computeTimestep(q, dx);

    computeFluxes(q, flux, dx, dt);

    // TODO - consider why these are the limits chosen
    for(int i = 1; i < cells + 1; ++i)
    {
      for(int var = 0; var < 3; ++var)
      {
	q[i][var] = q[i][var] - (dt/dx) * (flux[i+1][var] - flux[i][var]);
      }
    }

    t += dt;

    // TODO
    // Output some useful data to screen here so you know the code is running as expected
  }
 
  // Output

  std::ofstream rhoOutput("density.dat");
  std::ofstream velOutput("density.dat");
  std::ofstream preOutput("density.dat"); 

  for(unsigned int i=1 ; i < cells + 1 ; i++)
  {
    double x;
    // TODO
    // Compute the value of x associated with cell i
    
    std::array<double, 3> w = primitive(q[i]);
    
    rhoOutput << x << " " << w[0] << std::endl;
    velOutput << x << " " << w[1] << std::endl;
    preOutput << x << " " << w[2] << std::endl;
  }
  
  return 0;
}
