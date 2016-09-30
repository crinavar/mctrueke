//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  mctrueke                                                                    //
//  A multi-core implementation of the exchange Monte Carlo method.             //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  Copyright © 2015 Cristobal A. Navarro, Wei Huang.                           //
//                                                                              //
//  This file is part of mctrueke.                                              //
//  mctrueke is free software: you can redistribute it and/or modify            //
//  it under the terms of the GNU General Public License as published by        //
//  the Free Software Foundation, either version 3 of the License, or           //
//  (at your option) any later version.                                         //
//                                                                              //
//  mctrueke is distributed in the hope that it will be useful,                 //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of              //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               //
//  GNU General Public License for more details.                                //
//                                                                              //
//  You should have received a copy of the GNU General Public License           //
//  along with mctrueke.  If not, see <http://www.gnu.org/licenses/>.           //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////
#ifndef _STRUCTS_H_
#define _STRUCTS_H_

const char *symbols[] = {"E", "M", "sqE", "sqM", "quadM", "Xd", "F", "C", "X", "CORRLEN", "EXCHANGE", "ZSQE", "ZSQM"};
const char *filenames[] = {"energy.dat", "magnetization.dat", "sqenergy.dat", "sqmagnetization.dat", "quadmagnetization.dat", "dis_susceptibility.dat", "F.dat", "specific_heat.dat", "susceptibility.dat", "corrlen.dat", "exchange.dat", "zsqe.dat", "zsqm.dat"};

/* definitions for physical values */
#define NUM_PHYSICAL_VALUES 	13
#define NUM_SPECIAL				6

/* NORMAL: the following are computed at each ptstep */
#define E_POS 			0
#define M_POS 			1
#define SQE_POS			2
#define SQM_POS 		3
#define QUADM_POS 		4
#define Xd_POS 			5
#define F_POS		 	6
/* SPECIAL: the following are computed at each realization */
#define C_POS 			7
#define X_POS 			8
#define CORRLEN_POS		9
#define EXCHANGE_POS	10
#define ZSQE_POS		11
#define ZSQM_POS		12
#define C(x,y,z,L) 		((z)*(L)*(L)+(y)*(L)+(x))

/* forward struct declarations */
struct realization_data;
struct block_data;
struct mc_data;
struct obset;
struct setup;
struct float3;

typedef realization_data rdata_t;
typedef block_data bdata_t;
typedef mc_data mcdata_t;
typedef obset obset_t;
typedef setup setup_t;
typedef float3 float3_t;

// struct float3
struct float3{
    float x;
    float y;
    float z;
};

/* realization statistics data structure */
struct realization_data{
	int n;
	double mean;
	double stdev;
	double correlation;
	double avbstdev;
	double avbcorrelation;
	double w1;
	double w2;
	double x1;
	double lastx;	

};

/* block statistics data structure */
struct block_data{
	double mean;
	int n;
	double w1;
	double w2;
	double x1;
	double lastx;	
};

/* Monte Carlo step statistics data */
struct mc_data{
	double E;
	double sqE;
	double M;
	double sqM;
	double quadM;
	double F;
};

/* obervables set of values */
struct obset{
	rdata_t		rdata[NUM_PHYSICAL_VALUES];
	bdata_t 	bdata[NUM_PHYSICAL_VALUES];
	mcdata_t	mdata;
};

/* main structure */
struct setup{
	/* arguments / parameters */
	int L, N;
	float TR, dT, h;
	/* pt and simulation parameters */
	int pts, ds, ms, fs, cs, period;
	int blocks, realizations;
	int seed, oseed, R;
	int rthreads, sthreads;
    
    // distributed rng
    pcg32_random_t **rng;
    pcg32_random_t mainrng;
    
	int mzone;
	obset_t *obstable;
	/* temperature array, increasing */
	float *T;
	/* exchange array for each temp value */
	float *ex, *avex;
	/* index array for accessing a replica's temp and viceversa */
	int *rts, *trs;

	/* energy variables */
	float *exE, *E;
	/* magnetization variables */
	int *M;
	/* correlation length variables */
	float3_t *F1, *F2;

	/* set of host lattices and dH  */
	int **hlat, *hH;
	const char *obsfolder;
	const char *plotfolder;
};

#endif
