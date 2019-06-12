// FisherStefan.cpp
// Fisher-Stefan model solved by Raphson-Newton method
// Author: Maud El-Hachem
// June 2019
// QUT, Brisbane, Australia

#include <iomanip>  
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <cmath>        // std::abs

// table of the parameters needed for the simulation and used in the filename
// respectively in order:
// dx : Delta x
// dt : Delta t
// Tmax : time of simulation
// alpha: parameter alpha used in the initial conditions
// beta: parameter beta used in the initial conditions

// kappa: parameter in the Stefan condition
const double kappa = 0.45;

const double arg[7] = {1.0e-04,1.0e-03,20,0.5,1,kappa,1.0e-08};



// array of time steps for which to print the results
// for example, to print the results at time t=2, the time step is 2/dt
const int timetoPrint[] = {0*1000,4*1000,8*1000,12*1000,16*1000,20*1000}; 

const double pi = 3.1415926535897;
// this function implements Thomas algorithm to solve tridiagonal system Ax = d
// parameters: 	a, b, c	subdiagonal, diagonal and superdiagonal coefficients of the matrix A
// 				d		vector/array of right hand side values
//				x		solution vector
//				n		the size of the square matrix
void tridia(const double *a, const double *b, const double *c, const double *d, double *x, const int n)
{
	// allocate the memory
	double *bb, *dd, ff;
	bb = new double[n];
	dd = new double[n];

	for(int i=0; i<n; i++){
		bb[i] = b[i];
		dd[i] = d[i];
	}
	
	for(int i=1; i<n; i++){
		ff = a[i]/bb[i-1];
		bb[i] = bb[i]-c[i-1]*ff;
		dd[i] = dd[i]-dd[i-1]*ff;
	}

	int j = 0;
	x[n-1] = dd[n-1]/bb[n-1];
	for(int i = 0; i<n-1; i++){
		j = n-2-i;
		x[j] = (dd[j]-c[j]*x[j+1])/bb[j];
	}
	
	// free the memory
	delete [] bb;
	delete [] dd;
}

// this function prints the density solution, the current time, speed and boundary position
// parameters:  x			the current density
//				n			the number of nodes
//				T			the current time
//				c			the current wave speed
//				c			the current position of the moving boundary
void printSolutionTime(const double *x, const int n, const double T,const double c, const double L)
{
	std::ostringstream strs;
	for (int i = 0; i < 6; i++)
		strs << arg[i] << "_";
	strs << T << "_";
	std::string str = strs.str();
	std::string filename = str + "FisherStefan.bin";
	std::ofstream outfile (filename.c_str(), std::ofstream::binary);
	if (!outfile.is_open())
		return;
	outfile.write((char*)&T,sizeof(double));
	outfile.write((char*)&c,sizeof(double));
	outfile.write((char*)&L,sizeof(double));
  	for (int i = 0; i<n; i++)
		outfile.write((char*)&x[i],sizeof(double));
	outfile.close();
}

bool run()
{
	// VARIABLES DECLARATIONS
	bool success = false;

	// the coefficients arrays of the tridiagonal matrix
	// u_{i-1}^{j+1}, u_{i}^{j+1}, u_{i+1}^{j+1}, f^{j+1}
	double *ui_m, *ui, *ui_p, *f;
	// the density arrays u_i^{j+1}, u_i^{j} and Delta u
	double *u, *u_p, *delta;
	double dx = arg[0], dt = arg[1];
	double Tmax = arg[2];
	double alpha = arg[3], beta = arg[4];
	double tolerance = arg[6];

	double norme = 0;
	double L, L_p;
	
	// total number of nodes in the mesh
	int nodes = (int)round(beta/dx)+1;
	// total number of time steps for the simulation
	int steps = (int)round(Tmax/dt)+1;
	// travelling wave speed
	double c = 0;

	// allocate the memory
	ui_m = new double[nodes];
	ui = new double[nodes];
	ui_p = new double[nodes];
	f = new double[nodes];
	u = new double[nodes];  
	u_p = new double[nodes];  
	delta = new double[nodes]; 

	int p = 0; // iterator on the time in the array
	
	// initialise the density function with the initial conditions
	for(int i = 0; i<nodes-1; i++)
	{		
		u[i] = alpha;
		u_p[i] = alpha;
		delta[i] = 0;	
	}	
	u[nodes-1] = 0;
	u_p[nodes-1] = 0;
	delta[nodes-1] = 0;
	// initial moving boundary position
	L = beta; L_p = beta;	

	// step size in the fixed domain
	double dxi = dx/(beta);

	// NEWTON-RAPHSON METHOD
	// temporal integration loop
	for (int t=0; t<steps; t++){
		
		// print the current results
		if (t == timetoPrint[p])
		{
			printSolutionTime(u_p, nodes, t*dt, c, L_p);
		}

		// iterate until the norm <= tolerance
		do
		{	

			// left BC		
			ui_m[0] = 0;
			ui[0] = -1;
			ui_p[0] = 1;
			f[0] = -(u[1]-u[0]);
				
			//	right BC
			ui_m[nodes-1] = 0;
			ui[nodes-1] = 1;                 
			ui_p[nodes-1] = 0;
			f[nodes-1] = 0;
			
			// internal nodes
			for(int i = 1; i<nodes-1; i++)
			{	ui_m[i] = 1/(dxi*dxi*L*L) - i*dxi/L * (L-L_p)/(2*dt*dxi);
				ui[i] = -1/dt - 2/(dxi*dxi*L*L) + (1-2*u[i]);
				ui_p[i] = 1/(dxi*dxi*L*L) + i*dxi/L * (L-L_p)/(2*dt*dxi);
				f[i] = -(u[i-1] - 2*u[i] + u[i+1])/(dxi*dxi*L*L) 
						- i*dxi/L * (u[i+1] - u[i-1])* (L-L_p)/(2*dt*dxi) 
						- u[i]*(1-u[i]) + (u[i]-u_p[i])/dt;
			}

			// compute next time step values
			tridia(ui_m, ui, ui_p, f, delta, nodes);  	
			for(int i = 0; i<nodes; i++)
				u[i] += delta[i];
			L = L_p + dt*kappa*(u[nodes-2])/(L_p*dxi);
			
			// compute norm
			norme = 0;
			for(int i = 0; i<nodes; i++)
			{
				if (std::fabs(delta[i]) > norme)
					norme = std::fabs(delta[i]);
			}

		}while (norme > tolerance);
		
		// compute the speed for the current printing time
		if (t == timetoPrint[p]){
			c = (L-L_p)/dt;
			// go to the next printing time
			p++;
		}

		// update previous density and moving boundary position
		for(int i = 0; i<nodes; i++)
			u_p[i] = u[i];
		L_p = L;

	}

	delete [] ui_m;
	delete [] ui;
	delete [] ui_p;
	delete [] f;
	delete [] u;
	delete [] u_p;
	delete [] delta;
	return success;
}

int main(int argc, char **argv)
{
	bool success = run();
	return success;
}
