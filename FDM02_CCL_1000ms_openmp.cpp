/* Method 2 for 2D LR91,OpenMP support */
/* Operator splitting + FDM for PDE
   CCL for ODE */
/* ref: An Advanced Algorithm for Solving Partial Differential Equation in Cardiac Conduction. 1999. */
/* ref: Quadratic adaptive algorithm for solving cardiac action potential models. 2016*/
/* Xiang Zhou, 2017/11/30 */

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>  
#include <sys/timeb.h>

//*******FDM parameters for LR91 *******
int const nx = 1000, ny = 1000;//grid numbers
double dx = 0.015, dy = 0.015;//space step, 3cm*3cm
double D = 0.001;//D: diffusion coefficient cm^2/ms

/* Time Step */
double dt_max = 0.04; // Time step (ms)
double dt_min = 0.001;
double dt[nx + 1][nx + 1];
double t; // Time (ms)
int steps; // Number of Steps
int increment; // Loop Control Variable
int cutcount = 40 / dt_max;

int betmp = 0;
double dt_univ = dt_max;//ms universal maximum time-step siz
double D1V[nx + 2][nx + 2] = { 0 }, D2V[nx + 2][nx + 2] = { 0 }, dt_tmp = 0;
double Voffset = 0.1;//mv, a given membrane potential offset

/* Voltage */
double V[nx + 2][nx + 2]; // Initial Voltage (mv)
double dV2[nx + 2][nx + 2]; // second order derivatives of Voltage (mv)
double dVdt[nx + 2][nx + 2]; // first order derivatives of Voltage (mv)
double Vnew[nx + 2][nx + 2];// New Voltage (mV)
double dvdt; // Change in Voltage / Change in Time (mV/ms)
double dvdtnew; // New dv/dt (mV/ms)

/* Total Current and Stimulus */
double st; // Constant Stimulus (uA/cm^2)
double tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
int stimtime = (int)(0.6 / dt_max + 0.6); //Time period during which stimulus is applied 
double it[nx + 1][nx + 1]; // Total current (uA/cm^2)

/* Terms for Solution of Conductance and Reversal Potential */
const double R = 8314; // Universal Gas Constant (J/kmol*K)
const double frdy = 96485; // Faraday's Constant (C/mol)
double temp = 310; // Temperature (K)

/* Ion Concentrations */
double nai = 18;; // Intracellular Na Concentration (mM)
double nao = 140; // Extracellular Na Concentration (mM)
double cai[nx + 1][nx + 1]; // Intracellular Ca Concentration (mM)
double cao = 1.8; // Extracellular Ca Concentration (mM)
double ki = 145; // Intracellular K Concentration (mM)
double ko = 5.4; // Extracellular K Concentration (mM)

/* Fast Sodium Current (time dependant) */
double ina[nx + 1][ny + 1]; // Fast Na Current (uA/uF)
double gna = 23; // Max. Conductance of the Na Channel (mS/uF)
double ena = ((R*temp) / frdy)*log(nao / nai); // Reversal Potential of Na (mV)
//double am; // Na alpha-m rate constant (ms^-1)
//double bm; // Na beta-m rate constant (ms^-1)
//double ah; // Na alpha-h rate constant (ms^-1)
//double bh; // Na beta-h rate constant (ms^-1)
//double aj; // Na alpha-j rate constant (ms^-1)
//double bj; // Na beta-j rate constant (ms^-1)
//double mtau; // Na activation
//double htau; // Na inactivation
//double jtau; // Na inactivation
//double mss; // Na activation
//double hss; // Na inactivation
//double jss; // Na slow inactivation
double m[nx + 1][ny + 1]; // Na activation
double h[nx + 1][ny + 1]; // Na inactivation
double jj[nx + 1][ny + 1]; // Na slow inactivation

/* Current through L-type Ca Channel */
double isi[nx + 1][ny + 1]; // Slow inward current (uA/uF)
double esi[nx + 1][ny + 1]; // Reversal Potential of si (mV)
//double dcai; // Change in myoplasmic Ca concentration (mM)
//double ad; // Ca alpha-d rate constant (ms^-1)
//double bd; // Ca beta-d rate constant (ms^-1)
//double af; // Ca alpha-f rate constant (ms^-1)
//double bf; // Ca beta-f rate constant (ms^-1)

double d[nx + 1][ny + 1]; // Voltage dependant activation gate
double f[nx + 1][ny + 1]; // Voltage dependant inactivation gate
double fca[nx + 1][ny + 1]; // Ca dependant inactivation gate -from LR94
//double dss; // Steady-state value of activation gate d
//double taud; // Time constant of gate d (ms^-1)----mistake ???��ms��
//double fss; // Steady-state value of inactivation gate f
//double tauf; // Time constant of gate f (ms^-1)

/* Time-dependent potassium current*/
//double prnak = 0.01833;
//ek = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));
double gk = 0.282*sqrt(ko / 5.4); // Channel Conductance of Rapidly Activating K Current (mS/uF)
double ek = ((R*temp) / frdy)*log(ko / ki); // Reversal Potential of Rapidly Activating K Current (mV)
double ik[nx + 1][ny + 1]; // Rapidly Activating K Current (uA/uF)
double X[nx + 1][ny + 1]; // Rapidly Activating K time-dependant activation  --gate X in LR91
//double ax; // K alpha-x rate constant (ms^-1)
//double bx; // K beta-x rate constant (ms^-1)
//double xss; // Steady-state value of inactivation gate xr  --gate X in LR91
//double taux; // Time constant of gate xr (ms^-1) --gate X in LR91
//double Xi; // K time-independent inactivation --gate Xi in LR91

/* Potassium Current (time-independent) */
double ik1[nx + 1][ny + 1]; // Time-independent K current (uA/uF)
double gk1 = 0.6047*(sqrt(ko / 5.4)); // Channel Conductance of Time Independant K Current (mS/uF)
double ek1 = ((R*temp) / frdy)*log(ko / ki); // Reversal Potential of Time Independant K Current (mV)
//double ak1; // K alpha-ki rate constant (ms^-1)
//double bk1; // K beta-ki rate constant (ms^-1)
//double K1ss; // Steady-state value of K inactivation gate K1

/* Plateau Potassium Current */
double ikp[nx + 1][ny + 1]; // Plateau K current (uA/uF)
double gkp = 0.0183; // Channel Conductance of Plateau K Current (mS/uF)
double ekp = ek1; // Reversal Potential of Plateau K Current (mV)
//double kp; // K plateau factor

/* Background Current */
double ib[nx + 1][nx + 1]; // Background current (uA/uF)

//temporary gate value
double m0[nx + 1][nx + 1], h0[nx + 1][nx + 1], jj0[nx + 1][nx + 1];
double d0[nx + 1][nx + 1], f0[nx + 1][nx + 1];
double X0[nx + 1][nx + 1];

//performance compared
double Vmax, V_left = 0, V_right = 0, left_peak, right_peak, conduction_t = 0;
double APD90; // Time of 90% Repolarization 
double Vold, v_onset;

/* Ion Current Functions */
void comp_ina(int i, int j, double dt); // Calculates Fast Na Current
void comp_ical(int i, int j, double dt); // Calculates Currents through L-Type Ca Channel
void comp_ik(int i, int j, double dt); // Calculates Time-dependent K Current
void comp_ik1(int i, int j); // Calculates Time-Independent K Current
void comp_ikp(int i, int j); // Calculates Plateau K Current
void comp_ib(int i, int j); // Calculates Background Current
void comp_it(int i, int j); // Calculates Total Current
double get_it(int i, int j, double dt);
double CCL(int i, int j, double dt);
void comp_fluxes(int i, int j, double dt);
void Rush_Larsen(int i, int j, double dt);

FILE *single_ap;
void performance();

int main(int argc, char* argv[])
{
	/* Data File */
	FILE *ap;
	FILE *fevaluation;
	fevaluation = fopen("fevaluation", "w");
	single_ap = fopen("single_ap", "w");

	/* Time Loop Conditions */
	t = 0.0; // Time (ms)
	//	steps = (bcl*beats)/udt; // Number of ms
	st = -80.0; // Stimulus (mA)

	/* Beginning Ion Concentrations */
	nai = 18; // Initial Intracellular Na (mM)
	nao = 140; // Initial Extracellular Na (mM)
	ki = 145; // Initial Intracellular K (mM)
	ko = 5.4; // Initial Extracellular K (mM)
	//cai = 0.0002; // Initial Intracellular Ca (mM)
	cao = 1.8; // Initial Extracellular Ca (mM)

	int ncount, i, j;

	for (i = 1; i < nx + 1; i++){
		for (j = 1; j < ny + 1; j++){
			V[i][j] = -88.654973; // Initial Voltage (mv)
			m[i][j] = 0.000838;
			h[i][j] = 0.993336;
			jj[i][j] = 0.995484;
			d[i][j] = 0.000003;
			f[i][j] = 0.999745;
			X[i][j] = 0.000129;
			cai[i][j] = 0.0002; // Initial Intracellular Ca (mM)
			dt[i][j] = dt_max;
		}
	}

	int nstep = 10 / dt_max; // snapshot interval to save data files
	int index = 0;// filename index
	char filename[100];

	struct timeb start, end;
	int diff;
	ftime(&start);

	//	omp_set_num_threads(8);
	int tid = 0;
#pragma omp parallel private(ncount,tid)// create and destruction of threads takes time, so it's a good place to create all threads here!!
	{
		tid = omp_get_thread_num();// Get the thread ID
		printf("aaa from thread %d\n", tid);
		for (ncount = 1; ncount <= 1000 / dt_max; ncount++){//simulation time is 160ms
			if (tid == 0){// let thread 0 manage the boundaries
				for (i = 1; i < nx + 1; i++){
					//****no flux boundary conditions*****
					V[i][0] = V[i][1];
					V[i][ny + 1] = V[i][ny];
				}
				for (j = 1; j < ny + 1; j++){
					V[0][j] = V[1][j];
					V[nx + 1][j] = V[nx][j];
				}
				performance();
				//**** save data in file "ap"
				int fileflag = 0;

				for (i = 1; i < nx + 1; i++){
					for (j = 1; j < ny + 1; j++){
						if (ncount%nstep == 0){//
							if (fileflag == 0){
								sprintf(filename, "ap%d", index);
								ap = fopen(filename, "w");
								fileflag = 1;
								index++;
							}
							fprintf(ap, "%g\t", V[i][j]);
							if (j == ny){
								fprintf(ap, "\n");
							}
						}
					}
				}
				if (fileflag == 1){
					fclose(ap);
				}
				//**** save data in file "ap"
			}
#pragma omp barrier // wait for the boundary update

			//*********** step 1 *******
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					dV2[i][j] = D*((V[i + 1][j] + V[i - 1][j] - 2 * V[i][j]) / (dx*dx) + (V[i][j + 1] + V[i][j - 1] - 2 * V[i][j]) / (dy*dy));
				}
			}
#pragma omp barrier
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					//Forward Euler
					Vnew[i][j] = V[i][j] + dt_max / 2 * dV2[i][j];
					V[i][j] = Vnew[i][j];
				}
			}
#pragma omp barrier
			//*********** step 1 *******

			//*********** step 2 *******
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					it[i][j] = get_it(i, j, dt_max);
					dVdt[i][j] = -it[i][j];
				}
			}
#pragma omp barrier
			//*****stimulation with a plane waves****
			if (ncount >= 1 && ncount <= stimtime) { //stimulus is hold with 0.5 ms, 0.02*15 = 0.3 ms
#pragma omp for private(i,j) schedule(static)
				for (i = 1; i < nx + 1; i++){
					for (j = 1; j <= 5; j++){
						dVdt[i][j] = dVdt[i][j] + (-st);
					}
				}
#pragma omp barrier
			}

			if (ncount == 1){// in order to get D1V[i][j], for computing D2V[i][j] in CCL(i, j, dt_max);
				/* The first time step*/
#pragma omp for private(i,j) schedule(static)
				for (i = 1; i < nx + 1; i++){
					for (j = 1; j < ny + 1; j++){
						//Forward Euler
						D1V[i][j] = dVdt[i][j];
						comp_fluxes(i, j, dt[i][j]);
						Rush_Larsen(i, j, dt[i][j]);
						V[i][j] = V[i][j] + dt_max * dVdt[i][j];
					}
				}
#pragma omp barrier
			}
			else{
				double dt_sum;
#pragma omp for private(i,j, dt_sum) schedule(static)
				for (i = 1; i < nx + 1; i++){
					for (j = 1; j < ny + 1; j++){
						/******************* adjust or correct time step---CCL method ***********************/
						dt_sum = 0;
						dt[i][j] = CCL(i, j, dt_max);
						double dcai;
						do
						{
							dt_sum = dt_sum + dt[i][j];
							if (dt_sum<dt_max){
								D1V[i][j] = dVdt[i][j];
								comp_fluxes(i, j, dt[i][j]);
								Rush_Larsen(i, j, dt[i][j]);
								V[i][j] = V[i][j] + dt[i][j] * dVdt[i][j] + dt[i][j] * dt[i][j] * D2V[i][j] / 2;
							}
							else{
								dt[i][j] = dt_max - (dt_sum - dt[i][j]);// here is a new dt  !!!
								D1V[i][j] = dVdt[i][j];
								comp_fluxes(i, j, dt[i][j]);
								Rush_Larsen(i, j, dt[i][j]);
								V[i][j] = V[i][j] + dt[i][j] * dVdt[i][j] + dt[i][j] * dt[i][j] * D2V[i][j] / 2;
								break;
							}
							it[i][j] = get_it(i, j, dt[i][j]);
							if (ncount >= 1 && ncount <= stimtime && j >= 1 && j <= 5) {
								dVdt[i][j] = -it[i][j] + (-st);
							}
							else{
								dVdt[i][j] = -it[i][j];
							}
							dt[i][j] = CCL(i, j, dt[i][j]);
						} while (true);
					}
				}
				//!! using protective zone !!
				//if (t >= 10 + bcl*be && t <= (10 + 0.5 + bcl*be)){
				//	dt = 0.001;
				//}
				//dt_tmp = dt;
#pragma omp barrier
			}
			//*********** step 2 *******

			//*********** step 3 *******
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					dV2[i][j] = D*((V[i + 1][j] + V[i - 1][j] - 2 * V[i][j]) / (dx*dx) + (V[i][j + 1] + V[i][j - 1] - 2 * V[i][j]) / (dy*dy));
				}
			}
#pragma omp barrier
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					//Forward Euler
					Vnew[i][j] = V[i][j] + dt_max / 2 * dV2[i][j];
					V[i][j] = Vnew[i][j];
				}
			}
#pragma omp barrier
			//*********** step 3 *******

			if (tid == 0){
				t = t + dt_max;
			}

			//***********trancation 1/2 of the plane wave to generate a spiral wave******
			//if (ncount == cutcount){
			//	for (i = 1; i < nx / 2; i++){
			//		for (j = 1; j < ny; j++){
			//			V[i][j] = -88.654973; // Initial Voltage (mv)
			//			m[i][j] = 0.000838;
			//			h[i][j] = 0.993336;
			//			jj[i][j] = 0.995484;
			//			d[i][j] = 0.000003;
			//			f[i][j] = 0.999745;
			//			X[i][j] = 0.000129;
			//			cai[i][j] = 0.0002; // Initial Intracellular Ca (mM)
			//		}
			//	}
			//}
		}
	}
	ftime(&end);
	diff = (int)(1000.0*(end.time - start.time) + (end.millitm - start.millitm));
	conduction_t = (right_peak - left_peak)*0.001; //condution time from left side to right side
	fprintf(fevaluation, "%d\n%g\n%g\n", diff / 1000, APD90, dx*nx / conduction_t);
	fclose(fevaluation);
	fclose(single_ap);
}

//performance 
void performance(){
	fprintf(single_ap, "%g\n", V[nx / 2][ny / 2]);

	if (V[nx / 2][1] - V_left > 0){
		left_peak = t; // peak time at j=1
		V_left = V[nx / 2][1];
	}
	if (V[nx / 2][ny] - V_right > 0){
		right_peak = t; // peak time at j=ny
		V_right = V[nx / 2][ny];
	}
	if (V[nx / 2][ny / 2]>Vmax)
		Vmax = V[nx / 2][ny / 2];
	if (V[nx / 2][ny / 2] >= (Vmax - 0.9*(Vmax - (-88.654973))))
		APD90 = t; //  Time of 90% Repolarization 
}

double CCL(int i, int j, double dt){
	double dt_range;
	if (dt_univ > dt * 2){
		dt_range = dt * 2;
	}
	else{
		dt_range = dt_univ;
	}

	D2V[i][j] = (dVdt[i][j] - D1V[i][j]) / dt;
	double DiscriminantP = 0, DiscriminantN = 0, dtz = 0;
	if (dVdt[i][j] >= 0){
		DiscriminantP = dVdt[i][j] * dVdt[i][j] + 2 * D2V[i][j] * Voffset;
		if (D2V[i][j]>0){
			dt = (-dVdt[i][j] + sqrt(DiscriminantP)) / D2V[i][j];
		}
		else if (D2V[i][j]<0){
			dtz = -dVdt[i][j] / D2V[i][j];
			if (DiscriminantP >= 0){
				dt = (-dVdt[i][j] + sqrt(DiscriminantP)) / D2V[i][j];
			}
			else{
				dt = dtz;
			}
		}
	}
	else{
		DiscriminantN = dVdt[i][j] * dVdt[i][j] - 2 * D2V[i][j] * Voffset;
		if (D2V[i][j]>0){
			dtz = -dVdt[i][j] / D2V[i][j];
			if (DiscriminantN >= 0){
				dt = (-dVdt[i][j] - sqrt(DiscriminantN)) / D2V[i][j];
			}
			else{
				dt = dtz;
			}
		}
		else if (D2V[i][j]<0){
			dt = (-dVdt[i][j] - sqrt(DiscriminantN)) / D2V[i][j];
		}
	}

	if (dt>dt_range){
		dt = dt_range;
	}
	if (dt<dt_min){
		dt = dt_min;
	}
	return dt;
}

/********************************************************/
/* Functions that describe the currents begin here */

double get_it(int i, int j, double dt){
	comp_ina(i, j, dt);
	comp_ical(i, j, dt);
	comp_ik(i, j, dt);
	comp_ik1(i, j);
	comp_ikp(i, j);
	comp_ib(i, j);
	comp_it(i, j);
	return it[i][j];
}

/* Calculating the initial Ions Current*/
void comp_fluxes(int i, int j, double dt){
	double dcai = -0.0001*isi[i][j] + 0.07*(0.0001 - cai[i][j]);
	cai[i][j] = cai[i][j] + dcai*dt;//update Cai
}

/* update the gate value*/
void Rush_Larsen(int i, int j, double dt){
	//Fast sodium current
	/*gate variables can not be shared, should be local due to data racing !!!!!!!!*/
	double am = 0.32*(V[i][j] + 47.13) / (1 - exp(-0.1*(V[i][j] + 47.13)));
	double bm = 0.08*exp(-V[i][j] / 11);
	double ah, bh, aj, bj;
	if (V[i][j] < -40) {
		ah = 0.135*exp((80 + V[i][j]) / -6.8);
		bh = 3.56*exp(0.079*V[i][j]) + 310000 * exp(0.35*V[i][j]);
		aj = (-127140 * exp(0.2444*V[i][j]) - 0.00003474*exp(-0.04391*V[i][j]))*((V[i][j] + 37.78) / (1 + exp(0.311*(V[i][j] + 79.23))));
		bj = (0.1212*exp(-0.01052*V[i][j])) / (1 + exp(-0.1378*(V[i][j] + 40.14)));
	}
	else {
		ah = 0;
		bh = 1 / (0.13*(1 + exp((V[i][j] + 10.66) / -11.1)));
		aj = 0;
		bj = (0.3*exp(-0.0000002535*V[i][j])) / (1 + exp(-0.1*(V[i][j] + 32)));
	}
	double mtau = 1 / (am + bm);
	double htau = 1 / (ah + bh);
	double jtau = 1 / (aj + bj);

	double mss = am*mtau;
	double hss = ah*htau;
	double jss = aj*jtau;
	m[i][j] = mss - (mss - m[i][j])*exp(-dt / mtau);
	h[i][j] = hss - (hss - h[i][j])*exp(-dt / htau);
	jj[i][j] = jss - (jss - jj[i][j])*exp(-dt / jtau);


	//Slow inward current
	double ad = 0.095*exp(-0.01*(V[i][j] - 5)) / (1 + exp(-0.072*(V[i][j] - 5)));
	double bd = 0.07*exp(-0.017*(V[i][j] + 44)) / (1 + exp(0.05*(V[i][j] + 44)));
	double af = 0.012*exp(-0.008*(V[i][j] + 28)) / (1 + exp(0.15*(V[i][j] + 28)));
	double bf = 0.0065*exp(-0.02*(V[i][j] + 30)) / (1 + exp(-0.2*(V[i][j] + 30)));

	double taud = 1 / (ad + bd);
	double tauf = 1 / (af + bf);

	double dss = ad*taud;
	double fss = af*tauf;
	d[i][j] = dss - (dss - d[i][j])*exp(-dt / taud);
	f[i][j] = fss - (fss - f[i][j])*exp(-dt / tauf);

	//Time-dependent potassium current
	double ax = 0.0005*exp(0.083*(V[i][j] + 50)) / (1 + exp(0.057*(V[i][j] + 50)));
	double bx = 0.0013*exp(-0.06*(V[i][j] + 20)) / (1 + exp(-0.04*(V[i][j] + 20)));
	double taux = 1 / (ax + bx);
	double xss = ax*taux;
	X[i][j] = xss - (xss - X[i][j])*exp(-dt / taux);
}

//Fast sodium current
void comp_ina(int i, int j, double dt) {
	/*gate variables can not be shared, should be local due to data racing !!!!!!!!*/
	//double am = 0.32*(V[i][j] + 47.13) / (1 - exp(-0.1*(V[i][j] + 47.13)));
	//double bm = 0.08*exp(-V[i][j] / 11);
	//double ah, bh, aj, bj;
	//if (V[i][j] < -40) {
	//	ah = 0.135*exp((80 + V[i][j]) / -6.8);
	//	bh = 3.56*exp(0.079*V[i][j]) + 310000 * exp(0.35*V[i][j]);
	//	aj = (-127140 * exp(0.2444*V[i][j]) - 0.00003474*exp(-0.04391*V[i][j]))*((V[i][j] + 37.78) / (1 + exp(0.311*(V[i][j] + 79.23))));
	//	bj = (0.1212*exp(-0.01052*V[i][j])) / (1 + exp(-0.1378*(V[i][j] + 40.14)));
	//}
	//else {
	//	ah = 0;
	//	bh = 1 / (0.13*(1 + exp((V[i][j] + 10.66) / -11.1)));
	//	aj = 0;
	//	bj = (0.3*exp(-0.0000002535*V[i][j])) / (1 + exp(-0.1*(V[i][j] + 32)));
	//}
	//double mtau = 1 / (am + bm);
	//double htau = 1 / (ah + bh);
	//double jtau = 1 / (aj + bj);

	//double mss = am*mtau;
	//double hss = ah*htau;
	//double jss = aj*jtau;

	//m0[i][j] = mss - (mss - m[i][j])*exp(-dt / mtau);
	//h0[i][j] = hss - (hss - h[i][j])*exp(-dt / htau);
	//jj0[i][j] = jss - (jss - jj[i][j])*exp(-dt / jtau);

	ina[i][j] = gna*m[i][j] * m[i][j] * m[i][j] * h[i][j] * jj[i][j] * (V[i][j] - ena);
}

//Slow inward current
void comp_ical(int i, int j, double dt) {
	esi[i][j] = 7.7 - 13.0287*log(cai[i][j]);

	//double ad = 50 * 0.095*exp(-0.01*(V[i][j] - 5)) / (1 + exp(-0.072*(V[i][j] - 5)));
	//double bd = 50 * 0.07*exp(-0.017*(V[i][j] + 44)) / (1 + exp(0.05*(V[i][j] + 44)));
	//double af = 50 * 0.012*exp(-0.008*(V[i][j] + 28)) / (1 + exp(0.15*(V[i][j] + 28)));
	//double bf = 50 * 0.0065*exp(-0.02*(V[i][j] + 30)) / (1 + exp(-0.2*(V[i][j] + 30)));

	//double taud = 1 / (ad + bd);
	//double tauf = 1 / (af + bf);

	//double dss = ad*taud;
	//double fss = af*tauf;

	//d0[i][j] = dss - (dss - d[i][j])*exp(-dt / taud);
	//f0[i][j] = fss - (fss - f[i][j])*exp(-dt / tauf);

	isi[i][j] = 0.09*d[i][j] * f[i][j] * (V[i][j] - esi[i][j]);

	//dcai = -0.0001*isi[i][j] + 0.07*(0.0001 - cai[i][j]);

	//cai[i][j] = cai[i][j] + dcai*dt;
}

//Time-dependent potassium current
void comp_ik(int i, int j, double dt) {
	//double ax = 50 * 0.0005*exp(0.083*(V[i][j] + 50)) / (1 + exp(0.057*(V[i][j] + 50)));
	//double bx = 50 * 0.0013*exp(-0.06*(V[i][j] + 20)) / (1 + exp(-0.04*(V[i][j] + 20)));

	//double taux = 1 / (ax + bx);
	//double xss = ax*taux;

	//X0[i][j] = xss - (xss - X[i][j])*exp(-dt / taux);

	double Xi;
	if (V[i][j] > -100) {
		Xi = 2.837*(exp(0.04*(V[i][j] + 77)) - 1) / ((V[i][j] + 77)*exp(0.04*(V[i][j] + 35)));
	}
	else {
		Xi = 1;
	}

	ik[i][j] = gk*X[i][j] * Xi*(V[i][j] - ek);
}


//Time-independent potassium current
void comp_ik1(int i, int j) {
	double ak1 = 1.02 / (1 + exp(0.2385*(V[i][j] - ek1 - 59.215)));
	double bk1 = (0.49124*exp(0.08032*(V[i][j] - ek1 + 5.476)) + exp(0.06175*(V[i][j] - ek1 - 594.31))) / (1 + exp(-0.5143*(V[i][j] - ek1 + 4.753)));
	double K1ss = ak1 / (ak1 + bk1);

	ik1[i][j] = gk1*K1ss*(V[i][j] - ek1);
}

//Plateau potassium current
void comp_ikp(int i, int j) {
	double kp = 1 / (1 + exp((7.488 - V[i][j]) / 5.98));

	ikp[i][j] = gkp*kp*(V[i][j] - ekp);
}

//Background current
void comp_ib(int i, int j) {
	ib[i][j] = 0.03921*(V[i][j] + 59.87);
}

/* Total sum of currents is calculated here, if the time is between
stimtime = 0 and stimtime = 0.5 (ms), a stimulus is applied */
void comp_it(int i, int j) {
	//
	//	if (t >= 5 && t<(5 + 0.5)) {
	//		it[i][j] = st + ina[i][j] + isi[i][j] + ik[i][j] + ik1[i][j] + ikp[i][j] + ib[i][j];
	//	}else {
	it[i][j] = ina[i][j] + isi[i][j] + ik[i][j] + ik1[i][j] + ikp[i][j] + ib[i][j];
	//	}
}