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
#ifndef _TOOLS_H_
#define _TOOLS_H_

using namespace std;

/* print array function */
template <typename T>
void printarray(T *a, int n, const char *name);
/* forward function declarations */
void threadset(setup_t *s, int *tid, int *nt, int *r);
/* reset functions */
void reset_realization_statistics( setup_t *s, int tid, int a, int b);
void reset_mcmc_statistics( setup_t *s, int tid, int a, int b);
void reset_block_statistics( setup_t *s, int tid, int a, int b);
void reset(setup_t *s, int tid, int a, int b);
/* statistical averaging functions */
void accum_block_statistics( obset_t *obstable, int tid, int a, int b);
void accum_realization_statistics( setup_t *s, int a, int b, int realizations );
void accum_mcmc_statistics( setup_t *s, int tid, int a, int b); 

/* exchange temperatures function */
void extemp(setup_t *s, int r1, int r2);

// lattice up 
void set_lattice_up(int *s, int N){
for(int i=0; i < N; i++){
		s[i] = 1;
	}	
}

/* compare two floats */
int floatcomp(const void* elem1, const void* elem2){
    if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}

/* array reset */
template < typename T >
void reset_array(T *a, int n, T val){
	//printf("reseting array \n"); fflush(stdout);
	for(int i=0; i<n; ++i){
		a[i] = val;
	}
}
	
/* per realization reset */
void reset(setup_t *s, int tid, int a, int b){
	if(tid  == 0){
		printf("resets............."); fflush(stdout);
	}

#ifdef MEASURE
	/* reset block statistics */
	reset_block_statistics( s, tid, a, b);
#endif
	/* reset ex counters */
	reset_array<float>((float*)(s->ex + tid*(b-a)), b-a, 0.0f);

	/* reset average ex counters */
	reset_array<float>((float*)(s->avex + tid*(b-a)), b-a, 0.0f);

	/* reset index arrays */
	for(int i = a; i < b; ++i){
		s->rts[i] = s->trs[i] = i;
		s->avex[i] = 0.0;
        set_lattice_up(s->hlat[i], s->N);
	}
	if(tid  == 0){
		printf("ok\n"); fflush(stdout);
	}
}

/* set the thread indices for location */
void threadset(setup_t *s, int *tid, int *nt, int *r){
	/* get thread id */
	*tid = omp_get_thread_num();
	/* number of threads */
	*nt = omp_get_num_threads();
	/* 'r' replicas for each GPU */
	*r = (s->R)/(*nt);
	//printf("tid = %i   ngpus = %i\n   gpus[%i].i = %i\n", *tid, s->ngpus, s->gpus[*tid].i);
}

/* system call */
int run_sys_call(char *buffer){
    int res;
    res = system(buffer);
    if ( WEXITSTATUS(res) != 0 ) {
                syslog(LOG_CRIT," System call failed.\n");
                syslog(LOG_CRIT," %s\n",buffer);
    }
    return res;
}

/* make the output folders */
#ifdef MEASURE
void make_output_folders( const char *obs, const char *plot ){
	char command[256];
	sprintf(command, "mkdir -p %s/", obs);
	run_sys_call(command);
	sprintf(command, "mkdir -p %s/", plot);
	run_sys_call(command);
}
#endif

/* print average measures */
#ifdef MEASURE
void print_realization_statistics( setup_t *s ){
	for(int i = 0; i < s->R; ++i){
		printf("T[%i]=%.2f", i, s->T[i]);
		for(int j=0; j<NUM_PHYSICAL_VALUES; ++j){
			printf(" ,%s=[%.3f,%.3f,%.3f%%,%.2f]", symbols[j], s->obstable[i].rdata[j].mean, s->obstable[i].rdata[j].stdev, 100.0 * s->obstable[i].rdata[j].stdev / s->obstable[i].rdata[j].mean, s->obstable[i].rdata[j].correlation);
		}
		printf("\n");
	}
}
#endif

/* write the realization statistics */
#ifdef MEASURE
void write_realization_statistics( setup_t* s){
	FILE *fw;
	for(int j = 0; j < NUM_PHYSICAL_VALUES; j++){
		fw = fopen(string(string(s->obsfolder) + "/" + string(filenames[j])).c_str(), "w");
		//printf("plot file: %s\n", string(string(s->obsfolder) + "/" + string(filenames[j])).c_str());
		if(!fw){
				fprintf(stderr, "error opening file %s for writing\n", filenames[j]);
				exit(1);
		}
		fprintf(fw, "#T\t\t%s\t\tbstdev\t\tbsterr\t\tbcorr\t\trstdev\t\trsterr\t\trcorr\n", symbols[j]);
		/* temperature order */
		for(int i = 0; i < s->R; ++i){
			fprintf(fw, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", s->T[i], 
															s->obstable[i].rdata[j].mean, 
															s->obstable[i].rdata[j].avbstdev, 
															abs(s->obstable[i].rdata[j].avbstdev/sqrt(s->obstable[i].rdata[j].n)),
															abs(s->obstable[i].rdata[j].avbcorrelation), 
															s->obstable[i].rdata[j].stdev, 
															abs(s->obstable[i].rdata[j].stdev/sqrt(s->obstable[i].rdata[j].n)),
															s->obstable[i].rdata[j].correlation );

			fflush(fw);	
			//printf("writing replica %i\n", s->trs[i]);
		}
		fclose(fw);
	}
}
#endif

/* write the binder as : Qfrac = <m^4>/<m^2>^2 */
#ifdef MEASURE
void write_binder( setup_t *s ){
	FILE *fw;
	double x, y, dx, dy;
	fw = fopen(string(string(s->obsfolder) + "/" + string("binder.dat")).c_str(), "w");
	if(!fw){
		fprintf(stderr, "error opening file %s for writing\n", "binder.dat");
		exit(1);
	}
	fprintf(fw, "#T          B           dB\n");
	for(int r=0; r<s->R; ++r){
		x = s->obstable[r].rdata[QUADM_POS].mean;
		y = s->obstable[r].rdata[SQM_POS].mean;
		dx = s->obstable[r].rdata[QUADM_POS].stdev / sqrt((double)s->obstable[r].rdata[QUADM_POS].n);
		dy = s->obstable[r].rdata[SQM_POS].stdev / sqrt((double)s->obstable[r].rdata[SQM_POS].n);
		fprintf(fw, "%e\t%e\t%e\n", s->T[r], 	s->obstable[r].rdata[QUADM_POS].mean / (s->obstable[r].rdata[SQM_POS].mean * s->obstable[r].rdata[SQM_POS].mean),  (dx*y - dy*x) / (y*y));
	}
	fclose(fw);
}

void write_specific_heat( setup_t *s ){
	FILE *fw;
	double rsqT, dx, dy;
	fw = fopen(string(string(s->obsfolder) + "/" + string("Zspecific_heat.dat")).c_str(), "w");
	if(!fw){
		fprintf(stderr, "error opening file %s for writing\n", "Zspecific_heat.dat");
		exit(1);
	}
	fprintf(fw, "#T          C         dC\n");
	for(int r=0; r<s->R; ++r){
		rsqT = 1.0/((double)s->T[r] * (double)s->T[r]);
		dx = s->obstable[r].rdata[SQE_POS].stdev / sqrt((double)s->obstable[r].rdata[SQE_POS].n);
		dy = s->obstable[r].rdata[ZSQE_POS].stdev / sqrt((double)s->obstable[r].rdata[ZSQE_POS].n);
		fprintf(fw, "%e\t%e\t%e\n", s->T[r], 	(double)s->N * (s->obstable[r].rdata[SQE_POS].mean - s->obstable[r].rdata[ZSQE_POS].mean) * rsqT, abs(dx-dy));
	}
	fclose(fw);
}

void write_susceptibility( setup_t *s ){
	FILE *fw;
	double rT, dx, dy;
	fw = fopen(string(string(s->obsfolder) + "/" + string("Zsusceptibility.dat")).c_str(), "w");
	if(!fw){
		fprintf(stderr, "error opening file %s for writing\n", "Zsusceptibility.dat");
		exit(1);
	}
	fprintf(fw, "#T          X           dX\n");
	for(int r=0; r<s->R; ++r){
		rT = 1.0/((double)s->T[r]);
		dx = s->obstable[r].rdata[SQM_POS].stdev / sqrt((double)s->obstable[r].rdata[SQM_POS].n);
		dy = s->obstable[r].rdata[ZSQM_POS].stdev / sqrt((double)s->obstable[r].rdata[ZSQM_POS].n);
		fprintf(fw, "%e\t%e\t%e\n", s->T[r], 	(double)s->N * (s->obstable[r].rdata[SQM_POS].mean - s->obstable[r].rdata[ZSQM_POS].mean) * rT, abs(dx-dy));
	}
	fclose(fw);
}
#endif

/* print the magnetic field H */
void printH(int *h, int N){
	for(int i=0; i < N; i++)
		printf("%i \n", h[i]);
}

/* accumulate a variance step, for statistics */
void variance_step(double x,  int *n, double *mean, double *w1, double *w2, const double x1, double *lastx){
	
	*n = *n + 1;
	double d = x - *mean;
	double lastmean = *mean;

	*mean 	+= d/(*n);
	*w2 	+= d*(x-(*mean));

	// correlation
	*w1 += (x - *mean) * (*lastx - *mean) + (*n - 2) * (pow(*mean, 2.0) - pow(lastmean, 2.0)) + ((2.0*(*n-1)*(lastmean)) - x1 - *lastx) * (lastmean - *mean);
	*lastx = x;
}

/* reset per realization statistics */
void reset_realization_statistics( setup_t *s, int R){
	/* we do not need to use replica or temp order, it's just a reset */
	for(int r=0; r<R; ++r){
		for(int j=0; j<NUM_PHYSICAL_VALUES; j++){
			s->obstable[r].rdata[j].mean = 0.0;
			s->obstable[r].rdata[j].avbstdev = 0.0;
			s->obstable[r].rdata[j].avbcorrelation = 0.0;
			s->obstable[r].rdata[j].stdev = 0.0;
			s->obstable[r].rdata[j].correlation = 0.0;
			s->obstable[r].rdata[j].n = 0.0;
			s->obstable[r].rdata[j].w1 = 0.0;
			s->obstable[r].rdata[j].w2 = 0.0;
			s->obstable[r].rdata[j].x1 = 0.0;
			s->obstable[r].rdata[j].lastx = 0.0;
		}
	}
}

/* reset per block statistics */
void reset_block_statistics( setup_t *s, int tid, int a, int b){
	int q;
	for(int k = a; k < b; ++k){
		q = s->rts[k];
		for(int j=0; j<NUM_PHYSICAL_VALUES; j++){
			s->obstable[q].bdata[j].mean = 0.0;
			s->obstable[q].bdata[j].n = 0.0;
			s->obstable[q].bdata[j].w1 = 0.0;
			s->obstable[q].bdata[j].w2 = 0.0;
			s->obstable[q].bdata[j].x1 = 0.0;
			s->obstable[q].bdata[j].lastx = 0.0;
		}
	}
}

/* reset per mcmc statistics */
void reset_mcmc_statistics( setup_t *s, int tid, int a, int b){
	int q;
	for(int k = a; k < b; ++k){
		q = s->rts[k];
		s->obstable[q].mdata.E 		= 0.0;
		s->obstable[q].mdata.sqE 	= 0.0;
		s->obstable[q].mdata.M 		= 0.0;
		s->obstable[q].mdata.sqM 	= 0.0;
		s->obstable[q].mdata.quadM 	= 0.0;
		s->obstable[q].mdata.F		= 0.0;
	}
}

void metropolis(setup_t *s, int tid, int a, int b, const int ms, pcg32_random_t *trng){  
    for(int k = a; k < b; ++k){
        const double B = -2.0/s->T[s->trs[k]];
        #pragma omp parallel default(none) shared(s, tid, k) proc_bind(close) num_threads(s->sthreads) 
        {
            int ntid = omp_get_thread_num();
            for(int j=0; j < ms; ++j){
                //cpu_sweep(s->hlat[k], s->hH, s->h, s->L, B, tid, trng);
                multicore_sweep(s->hlat[k], s->hH, s->h, s->L, B, tid, &(s->rng[tid][ntid]), s->sthreads);
            }
        }
    }
}

/* measured mcmc simulation */
void simulation(setup_t *s, int tid, int a, int b, pcg32_random_t *trng){
	for(int i = 0; i < s->fs; i++){
		/* use replica order */
		for(int k = a; k < b; ++k){
			cpu_sweep(s->hlat[k], s->hH, s->h, s->L, -2.0/s->T[s->trs[k]], tid, trng);
		}
#ifdef MEASURE
		/* accumulate mcmc statistics */
		accum_mcmc_statistics(s, tid, a, b);
#endif
	}
}

/* accum per block measures. this code runs inside an openmp parallel pragma */
#ifdef MEASURE
void accum_mcmc_statistics( setup_t *s, int tid, int a, int b){
	double E, M, F;
	double invN = 1.0/(double)s->N;
	int L = s->L;
	/* accumulate (write) data in temp order, but read from replica order */
	for(int k = a; k < b; ++k){
		int q = s->trs[k];

		/* CPU measurements with double and see any real improvements for using double precision */
		E 	= (double)computeE(s->hlat[k], s->hH, s->h, L)*invN;
		M 	= abs((double)computeM(s->hlat[k], s->N))*invN;
		F 	= (double)computeF(s->hlat[k], L, L, L)*invN;

		/* write in temp order */
		s->obstable[q].mdata.E		+= E;
		s->obstable[q].mdata.sqE 	+= E*E;
		s->obstable[q].mdata.M		+= M;
		s->obstable[q].mdata.sqM 	+= M*M;
		s->obstable[q].mdata.quadM	+= M*M*M*M;
		s->obstable[q].mdata.F		+= F;
	}
}
#endif

/* accumulate per block statistics */
#ifdef MEASURE
void accum_block_statistics( setup_t *s, int tid, int a, int b ){
	int q;
	double invsteps = 1.0/(double)s->fs;
	/* array of observables, ommiting the special ones */
	double values[NUM_PHYSICAL_VALUES - NUM_SPECIAL];
	double N = (double)s->N;
	/* use temp order */
	for(int k = a; k < b; ++k){
		q = s->trs[k];
		values[E_POS] 	= s->obstable[q].mdata.E * invsteps;
		values[M_POS] 	= s->obstable[q].mdata.M * invsteps;
		values[SQE_POS] = s->obstable[q].mdata.sqE * invsteps;
		values[SQM_POS] = s->obstable[q].mdata.sqM * invsteps;
		values[QUADM_POS] = s->obstable[q].mdata.quadM * invsteps;
		values[Xd_POS] 	= N * (values[M_POS] * values[M_POS]) * invsteps;
		values[F_POS] 	= s->obstable[q].mdata.F * invsteps;
		for(int j=0; j<NUM_PHYSICAL_VALUES - NUM_SPECIAL; j++){
			if( s->obstable[q].bdata[j].n == 0 ){
				s->obstable[q].bdata[j].x1 = values[j];
			}
			variance_step(	values[j], 	&(s->obstable[q].bdata[j].n), 		&(s->obstable[q].bdata[j].mean),	&(s->obstable[q].bdata[j].w1), &(s->obstable[q].bdata[j].w2), 	s->obstable[q].bdata[j].x1, 		&(s->obstable[q].bdata[j].lastx) );
		}
	}
	//for(int k = a; k < b; ++k){
	//	printf("R%i  MEAN BLOCK        E = %f       M = %f\n\n", k, s->obstable[k].bdata[E_POS].mean, s->obstable[k].bdata[M_POS].mean);
	//}
}
#endif

/* accum realization statistics */
#ifdef MEASURE
void accum_realization_statistics( setup_t *s, int tid, int a, int b, int realizations ){
	for(int k = a; k < b; ++k){
		/* traverse replicas, but access in temp order */
		int q = s->trs[k];
		/* the special observables are ommited here */
		for(int j=0; j<NUM_PHYSICAL_VALUES - NUM_SPECIAL; ++j){
			if( s->obstable[q].rdata[j].n == 0 ){
				s->obstable[q].rdata[j].x1 = s->obstable[q].bdata[j].mean;
			}
			variance_step(s->obstable[q].bdata[j].mean, &(s->obstable[q].rdata[j].n), &(s->obstable[q].rdata[j].mean),&(s->obstable[q].rdata[j].w1), &(s->obstable[q].rdata[j].w2),
							s->obstable[q].rdata[j].x1,&(s->obstable[q].rdata[j].lastx) );

			s->obstable[q].rdata[j].avbstdev 				+= sqrt( s->obstable[q].bdata[j].w2 / ((double)s->obstable[q].bdata[j].n - 1.0)) / (double)realizations;
			s->obstable[q].rdata[j].avbcorrelation 		+= (s->obstable[q].bdata[j].w1 / s->obstable[q].bdata[j].w2) / (double)realizations;
		}
		/* auxiliary variables for computing observables */
		double A, sqA, N, L, val, F, T;
		L = (double)s->L;
		N = (double)s->N;
		T = (double)s->T[q];



		/* specific heat */
		A = s->obstable[q].bdata[E_POS].mean;
		sqA = s->obstable[q].bdata[SQE_POS].mean;
		val = N * (sqA - A*A) / (T*T); 
		if( s->obstable[q].rdata[C_POS].n == 0 ){
			s->obstable[q].rdata[C_POS].x1 = val;
		}
		variance_step(val, &(s->obstable[q].rdata[C_POS].n), &(s->obstable[q].rdata[C_POS].mean), &(s->obstable[q].rdata[C_POS].w1), &(s->obstable[q].rdata[C_POS].w2), s->obstable[q].rdata[C_POS].x1, 
						&(s->obstable[q].rdata[C_POS].lastx));

		
		/* susceptibility */
		A = s->obstable[q].bdata[M_POS].mean;
		sqA = s->obstable[q].bdata[SQM_POS].mean;
		val = N * (sqA - A*A) / T; 
		if( s->obstable[q].rdata[X_POS].n == 0 ){
			s->obstable[q].rdata[X_POS].x1 = val;
		}
		variance_step(val, &(s->obstable[q].rdata[X_POS].n), &(s->obstable[q].rdata[X_POS].mean), &(s->obstable[q].rdata[X_POS].w1), &(s->obstable[q].rdata[X_POS].w2), s->obstable[q].rdata[X_POS].x1, 
						&(s->obstable[q].rdata[X_POS].lastx));

		/* corrlen */
		sqA = s->obstable[q].bdata[SQM_POS].mean;
		F = s->obstable[q].bdata[F_POS].mean;
		//printf("sqM = %f        F = %f          sqM/F = %f\n", sqA, F, sqA/F);
		if( sqA >= F ){
			val = sqrt( sqA/F - 1.0 ) / L;
			if( s->obstable[q].rdata[CORRLEN_POS].n == 0 ){
				s->obstable[q].rdata[CORRLEN_POS].x1 = val;
			}
			variance_step(val, &(s->obstable[q].rdata[CORRLEN_POS].n), &(s->obstable[q].rdata[CORRLEN_POS].mean), &(s->obstable[q].rdata[CORRLEN_POS].w1), &(s->obstable[q].rdata[CORRLEN_POS].w2),	
							s->obstable[q].rdata[CORRLEN_POS].x1, &(s->obstable[q].rdata[CORRLEN_POS].lastx));
		}					

		/* ZSQE */
		val = pow(s->obstable[q].bdata[E_POS].mean, 2.0);
		if( s->obstable[q].rdata[ZSQE_POS].n == 0 ){
			s->obstable[q].rdata[ZSQE_POS].x1 = val;
		}
		variance_step(val, &(s->obstable[q].rdata[ZSQE_POS].n), &(s->obstable[q].rdata[ZSQE_POS].mean), &(s->obstable[q].rdata[ZSQE_POS].w1), &(s->obstable[q].rdata[ZSQE_POS].w2),	
						s->obstable[q].rdata[ZSQE_POS].x1, &(s->obstable[q].rdata[ZSQE_POS].lastx));

		/* ZSQM */
		val = pow(s->obstable[q].bdata[M_POS].mean, 2.0);
		if( s->obstable[q].rdata[ZSQM_POS].n == 0 ){
			s->obstable[q].rdata[ZSQM_POS].x1 = val;
		}
		variance_step(val, &(s->obstable[q].rdata[ZSQM_POS].n), &(s->obstable[q].rdata[ZSQM_POS].mean), &(s->obstable[q].rdata[ZSQM_POS].w1), &(s->obstable[q].rdata[ZSQM_POS].w2),	
						s->obstable[q].rdata[ZSQM_POS].x1, &(s->obstable[q].rdata[ZSQM_POS].lastx));

		/* exchange rates */
		val = 2.0 * s->avex[k];
		if(s->obstable[q].rdata[EXCHANGE_POS].n == 0){
			s->obstable[q].rdata[EXCHANGE_POS].x1 = val;
		}
		variance_step(val, &(s->obstable[q].rdata[EXCHANGE_POS].n), &(s->obstable[q].rdata[EXCHANGE_POS].mean), &(s->obstable[q].rdata[EXCHANGE_POS].w1), &(s->obstable[q].rdata[EXCHANGE_POS].w2),
							s->obstable[q].rdata[EXCHANGE_POS].x1, &(s->obstable[q].rdata[EXCHANGE_POS].lastx));
	}
}
#endif

/* produce realization statistics */
void make_realization_statistics( setup_t *s ){
	/* use replica order */
	for(int i=0; i<s->R; ++i){
		for(int j=0; j<NUM_PHYSICAL_VALUES; j++){
			s->obstable[i].rdata[j].stdev 			= sqrt(s->obstable[i].rdata[j].w2 / ((double)s->obstable[i].rdata[j].n - 1.0));
			s->obstable[i].rdata[j].correlation 	= abs((s->obstable[i].rdata[j].w1 / s->obstable[i].rdata[j].w2));
		}
	}
}

/* CPU random configuration for the lattice */
void random_configuration(setup_t *s, int N, int* lat, pcg32_random_t *trng){
	for(int i=0; i < N; i++){
		if(pcg32randn(trng) >= 0.5)
			lat[i] =  1;
		else
			lat[i] = -1;
	}	
}

/* CPU one configuration for the lattice */
void one_configuration(int N, int* lat){
	for(int i=0; i < N; i++){
			lat[i] =  1;
	}	
}

template <typename T>
void printarray(T *a, int n, const char *name){
	cout << name << "\t = [";
	for(int i = 0; i < n; ++i){
		cout << a[i] << ", ";
	}
	printf("]\n");
}

template <typename T>
void printindexarray(T *a, int *ind, int n, const char *name){
	cout << name << "\t = [";
    for(int i = 0; i < n; ++i){
        cout << a[ind[i]] << ", ";
    }
	printf("]\n");
}

/* CPU random magnetic field H */
void hdist(int *hH, const int N, const int tid, const int a, const int b, pcg32_random_t *trng){
	/* generate dist using CPU */
    if(tid == 0){
        random_Hi(hH, N, trng);
    }
    #pragma omp barrier
}

/* physical results */
#ifdef MEASURE
void physical_results(setup_t *s){
	make_realization_statistics(s);
	/* write special average measures */
	write_binder(s);
	write_specific_heat(s);
	write_susceptibility(s);
	/* write average measures */
	write_realization_statistics(s);
}
#endif

/* free memory for CPU and GPU */
void freemem(setup_t *s){
	for(int i = 0; i < s->R; ++i){
		//printf("feeing replica data %i.....", r); fflush(stdout);
		free(s->hlat[i]);
		//printf("ok\n"); fflush(stdout);
	}
	free(s->hlat);
	free(s->E);
	free(s->exE);
	free(s->hH);
	free(s->M);
	free(s->F1);
	free(s->F2);
	free(s->ex);
	free(s->avex);
#ifdef MEASURE
	free(s->obstable);
#endif
}

#endif
