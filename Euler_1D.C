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

// Calculate energy using equations on slide 14 of Euler lecture notes
double calc_E(double r, double u, double p, double gamma) {
    // where r is density, u is velocity, and p is pressure of an ideal gas
    return p/(gamma-1) + (0.5 * r * u * u);
}


std::array<double, 3> array_addition(std::array<double, 3> firstArray, std::array<double, 3> secondArray)
{

  std::array<double, firstArray.size() > outputArray;
  for (int i = 0; i < firstArray.size(); i++)
  {
    outputArray[i] = firstArray[i] + secondArray[i];
  }
  return outputArray;
}

std::array<double, 3> array_subtraction(std::array<double, 3> firstArray, std::array<double, 3> secondArray)
{

  std::array<double, firstArray.size() > outputArray;
  for (int i = 0; i < firstArray.size(); i++)
  {
    outputArray[i] = firstArray[i] - secondArray[i];
  }
  return outputArray;
}

// Convert from conserved variables to primitive variables
std::array<double, 3> primitive(std::array<double, 3>& q_i, double gamma){
  
  std::array<double, 3> w_i;

  // Convert the conserved state variables q_i to the primitive variables w_i
  w_i[0] = q_i[0];           // density
  w_i[1] = q_i[1] / q_i[0];  // velocity
  w_i[2] = (gamma-1) * (q_i[2] - 0.5 * q_i[0] * w_i[1] * w_i[1]);  // pressure
  
  return w_i;
}

// Convert primitive state vector to conservative state vector
std::array<double, 3> conservative(std::array<double, 3>& w_i, double gamma)
{
  std::array<double, 3> q_i;

  // Convert the primitive variables w_i to the conserved variables q_i
  // w_i = (r, u, p) and q_i = (r, r * u, E) - from slides 11 and 18

  q_i[0] = w_i[0];                                 // density
  q_i[1] = w_i[0] * w_i[1];                        // momentum
  q_i[2] = calc_E(w_i[0], w_i[1], w_i[2], gamma);  // energy

  return q_i;
}

// Fill the conservative state vector with the initial data
void initialiseData(double gamma, std::vector<std::array<double, 3> >& q, int test, double& finalT)
{
  int n = q.size();

  if(test == 1)
  {
    finalT = 0.25;

    for(int i = 0; i < n; ++i)
    {
      std::array<double, 3> w;
      if(i < n/2)
      {
	w = {1.0,0.0,1.0};
      }
      else
      {
	w = {0.125,0.0,0.1};
      }
      
      std::array<double, 3> q_i = conservative(w, gamma);
      q.at(i)[0] = q_i[0];
      q.at(i)[1] = q_i[1];
      q.at(i)[2] = q_i[2];
    }
  }
  else if(test == 2)
  {
    finalT = 0.25;

    for(int i = 0; i < n; ++i)
    {
      std::array<double, 3> w;
      if(i < n/2)
      {
  w = {1.0,0.0,1.0};
      }
      else
      {
  w = {0.125,0.0,0.1};
      }
      
      std::array<double, 3> q_i = conservative(w, gamma);
      q[i][0] = q_i[0];
      q[i][1] = q_i[1];
      q[i][2] = q_i[2];
    }
  }
  else if(test == 3)
  {
    finalT = 0.25;

    for(int i = 0; i < n; ++i)
    {
      std::array<double, 3> w;
      if(i < n/2)
      {
  w = {1.0,0.0,1.0};
      }
      else
      {
  w = {0.125,0.0,0.1};
      }
      
      std::array<double, 3> q_i = conservative(w, gamma);
      q[i][0] = q_i[0];
      q[i][1] = q_i[1];
      q[i][2] = q_i[2];
    }
  }
  else if(test == 4)
  {
    finalT = 0.25;

    for(int i = 0; i < n; ++i)
    {
      std::array<double, 3> w;
      if(i < n/2)
      {
  w = {1.0,0.0,1.0};
      }
      else
      {
  w = {0.125,0.0,0.1};
      }
      
      std::array<double, 3> q_i = conservative(w, gamma);
      q[i][0] = q_i[0];
      q[i][1] = q_i[1];
      q[i][2] = q_i[2];
    }
  }
  else if(test == 5)
  {
    finalT = 0.25;

    for(int i = 0; i < n; ++i)
    {
      std::array<double, 3> w;
      if(i < n/2)
      {
  w = {1.0,0.0,1.0};
      }
      else
      {
  w = {0.125,0.0,0.1};
      }
      
      std::array<double, 3> q_i = conservative(w, gamma);
      q[i][0] = q_i[0];
      q[i][1] = q_i[1];
      q[i][2] = q_i[2];
    }
  }
  else
  {
    // TODO
    // Add suitable error handling
  }
}

// Compute flux-vector corresponding to given state vector
// Equations for flux given on slide 11 of Euler notes
std::array<double, 3> flux(std::array<double, 3>& q_i, double gamma)
{
  std::array<double, 3> f;
  std::array<double, 3> w_i = primitive(q_i, gamma);

  f[0] = w_i[0] * w_i[1];                    // r * u
  f[1] = w_i[0] * w_i[1] * w_i[1] + w_i[2];  // r * u^2 + p
  f[2] = (q_i[2] + w_i[2]) * w_i[1];         // (E + p) * u
  
  return f;
}

// Compute the maximum stable time-step for the given data
double computeTimestep(double gamma, int cellNumber, std::vector<std::array<double, 3>> &q, double dx, double CFL)
{

  double dt = 0, alphaOld = 0, alphaMax = 0;
  std::vector<double> Cs;
  std::cout << "Gamma: " << gamma << std::endl;
  std::cout << "cellNumber: " << cellNumber << std::endl;
  std::cout << "dx: " << dx << std::endl;
  std::cout << "CFL: " << CFL << std::endl;

  // Compute the maximum wavespeed over the entire domain, and use this to compute the timestep
  for(int i=0 ; i < cellNumber ; i++)
  {
  /* w_i = (r, u, p) */
    Cs.push_back(sqrt( gamma * (q[2][i]/q[0][i]) ));
    std::cout << "Cs: " << Cs[i] << std::endl;
    alphaOld = abs(q.at(1)[i]) + Cs[i];
    std::cout << "alpha: " << alphaOld << std::endl;

    if (alphaOld > alphaMax) alphaMax = alphaOld;
  }

  dt = (CFL*dx)/alphaMax;
  std::cout << "dt: " << dt << std::endl;

  return dt;
}

// Compute the FORCE flux between two states uL and uR in coordinate direction coord.
std::array<double, 3> FORCEflux(double gamma, std::array<double, 3>& q_iMinus1, std::array<double, 3>& q_i, double dx, double dt)
{

  std::array<double, 3> fluxLF_temp1;
  std::array<double, 3> fluxLF_temp2;

  // Compute the Richtmyer flux using q_i and q_iMinus1
  std::array<double, 3> fluxRM;
  std::array<double, 3> q_iHalf;
  double halfDeltaTDeltaX = 0.5 * (dt/dx);

  fluxLF_temp1 = array_addition(q_iMinus1, q_i);
  fluxLF_temp2 = array_subtraction(flux(q_iMinus1, gamma), flux(q_i, gamma));

  // q_iHalf = 0.5 * (array_addition(q_iMinus1, q_i)) + halfdelta * 
  for (int i=0;i < fluxLF_temp1.size(); i++){
    fluxLF_temp1[i] *= 0.5;
    fluxLF_temp2[i] *= halfDeltaTDeltaX;
  }
  q_iHalf = array_addition(fluxLF_temp1, fluxLF_temp2);
  fluxRM = flux(q_iHalf, gamma);

  // Compute the Lax-Friedrichs flux using q_i and q_iMinus1
  std::array<double, 3> fluxLF;

  double halfDeltaXDeltaT = 0.5 * (dx/dt);
  // fluxLF = flux_temp1 + flux_temp2;
  fluxLF_temp1 = array_addition(flux(q_iMinus1, gamma), flux(q_i, gamma));
  fluxLF_temp2 = array_subtraction(q_iMinus1, q_i);

  for (int i=0;i < fluxLF_temp1.size(); i++){
    fluxLF_temp1[i] *= 0.5;
    fluxLF_temp2[i] *= halfDeltaXDeltaT;
  }

  fluxLF = array_addition(fluxLF_temp1, fluxLF_temp2);

  std::array<double, 3> fluxForce;
  // Compute the FORCE flux
  // fluxForce =  0.5 * (fluxRM + fluxLF);

  for (int i=0; i<3; i++){
    fluxForce[i] = 0.5 * (fluxLF[i] + fluxRM[i]);
  }
  
  return fluxForce;
}

// Compute the array of fluxes from the given data array
void computeFluxes(double gamma, std::vector<std::array<double, 3> >& q, std::vector<std::array<double, 3> >& flux, double dx, double dt)
{
  int n = q.size();
  std::cout << "N: " << n << std::endl;
  // Flux for current cell calculated using values from previous cell. Don't have values for cell i=-1 so can't start loop at i=0
  for(unsigned int i=1 ; i < n ; i++)
  {
    flux[i] = FORCEflux(gamma, q[i-1], q[i], dx, dt);
  }
}

int main(void)
{
  // User-defined parameters
  int cells;
  int test;

  std::cout << "Number of cells: " << std::flush;
  std::cin >> cells;

  std::cout << "Test number (1-5): " << std::flush;
  std::cin >> test;
  

  const double gamma = 1.4;
  const double CFL = 0.9;
  double finalT;
  
  // Set initial vectors:
  std::vector<std::array<double, 3> > q(cells+2);
  std::vector<std::array<double, 3> > flux(cells+2);
  // TODO - consider the extent of the vectors, the flux vector actually contains one more cell than it needs to
  
  initialiseData(gamma, q, test, finalT);
  
  // Do simulation
  double t = 0;
  double dx = 1.0 / cells;
  // You may wish to change the 1.0 to a specification of your own domain size

  std::cout << "Simulation started..." << std::endl;
  
  while( t < finalT )
  {
    // TODO
    // Implement this function, and then uncomment
    // computeDomainBoundaries(q)
    
    double dt = computeTimestep(gamma, cells, q, dx, CFL);

    computeFluxes(gamma, q, flux, dx, dt);

    // TODO - consider why these are the limits chosen
    for(int i = 1; i < cells + 1; ++i)
    {
      for(int var = 0; var < 3; ++var)
      {
        q[i][var] = q[i][var] - (dt/dx) * (flux[i+1][var] - flux[i][var]);
      } 
    }

    t += dt;

    // Output some useful data to screen here so you know the code is running as expected
    std::cout << "    t = " << t << " => " << t/finalT * 100.0 << "%" << std::endl;
  }
 
  // Output

  std::ofstream rhoOutput("density.dat");
  std::ofstream velOutput("velocity.dat");
  std::ofstream preOutput("pressure.dat"); 

  for(unsigned int i=1 ; i < cells + 1 ; i++) {

    double x = (double(i) + 0.5) * dx;
    std::array<double, 3> w = primitive(q[i], gamma);
    
    rhoOutput << x << " " << w[0] << std::endl;
    velOutput << x << " " << w[1] << std::endl;
    preOutput << x << " " << w[2] << std::endl;

  }
  
  return 0;
}