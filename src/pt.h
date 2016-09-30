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


#ifndef _PT_H_
#define _PT_H_

/* forward function declarations */
void ptenergies(setup_t *s, int tid, int a, int b);
int exchange(setup_t *s, int tid, int a, int b, int p, pcg32_random_t *trng);
void measure(setup_t *s, int tid, int a, int b, int p, pcg32_random_t *trng);
void pt(setup_t *s, int tid, int a, int b, pcg32_random_t *trng);
void swap(setup_t *s, int a, int b );

/* pt(): parallel tempering main loop */
void pt(setup_t *s, int tid, int a, int b, pcg32_random_t *trng){
        double t1;
		/* reset ex counters */
		reset_array<float>((float*)(s->ex + tid*(b-a)), b-a, 0.0f);
		/* reset average ex counters */
		reset_array<float>((float*)(s->avex + tid*(b-a)), b-a, 0.0f);
		/* average exchanges */
		double avex = 0.0;
		/* progress printing */
		if( tid == 0 ){
			printf("ptsteps............0%%"); fflush(stdout);
            t1 = omp_get_wtime();
		}
		/* parallel tempering */
		for(int p = 0; p < s->pts; ++p){
			/* pt metropolis phase */
			metropolis(s, tid, a, b, s->ms, trng);
			/* compute energies for exchange */
			ptenergies(s, tid, a, b);
			#ifdef MEASURE
				/* measure */
				if((p % s->period) == 0){
					//printf("measuring at ptstep = %i    period = %i\n", p, s->period);
					measure(s, tid, a, b, p, trng);
				}
			#endif
			/* exchange phase */
			avex += (double)exchange(s, tid, a, b, p, trng)/(double)s->pts;
			/* progress printing */
			if(tid == 0){
				printf("\rptsteps............%i%%", (100*(p+1))/s->pts); fflush(stdout);
				//printf("\navex = %f\n", avex);
			}
		}
		/* compute the average exchange for each replica r_i and r_i-1 */
		for(int i = a; i < b; ++i){
			s->avex[i] = s->ex[i] / (double)s->pts;
		}
		/* progress printing */
		if( tid  == 0 ){
			//printf(" %.3fs ", sdkGetTimerValue(&(s->timer))/1000.0f);
			printf("\t[<ex> = %.3f] %.2f secs\n\n", avex / ((double)(s->R-1)/2.0), omp_get_wtime() - t1);
            //printindexarray<float>(s->exE, s->rts, s->R, "exE");
            //printarray<float>(s->avex, s->R, "avex");
		}
}

/* energies for all replicas using Kepler shfl reductions and streams */
void ptenergies(setup_t *s, int tid, int a, int b){
	/* compute one energy reduction for each replica */
	for(int k = a; k < b; ++k){
		s->exE[k]   = (double)computeE(s->hlat[k], s->hH, s->h, s->L);
		//printf("        cpu E[%i] = %f\n", k, s->exE[k]);
	}
    #pragma omp barrier
	//if(tid == 0){
    //    printarray<float>(s->exE, s->R, "exE");
    //}
    //getchar();
	//printf("a = %i    b = %i\n", a, b);
	//for(int k = a; k < b; ++k){
	//	printf("compute E[%i] = %f\n", k, s->E[k]);
	//}
}

/* swap temperatures */
void swap(setup_t *s, int a, int b ){
	//printf("\nswapping T%i  T%i\n", a, b);
	int t1, t2, taux, raux;
	t1 = s->rts[a];
	t2 = s->rts[b];
	taux = s->trs[t1];
	raux = s->rts[a];

	/* swap rts */
	s->rts[a] = s->rts[b];
	s->rts[b] = raux; 

	/* swap trs */
	s->trs[t1] = s->trs[t2];
	s->trs[t2] = taux;
}

/* exchange phase */
int exchange(setup_t *s, int tid, int a, int b, int p, pcg32_random_t *trng){
	/* count the number of exchanges */
	int ex = 0;
	/* sync openmp threads before entering the exchange phase */
	#pragma omp barrier
	if(tid == 0){
		double delta = 0.0;
		/* traverse in reverse temperature order */
		//printf("\n");
		for(int k = s->R-1; k > 0; --k){						
				/* alternate between odd and even replicas */
				if((k % 2) == (p % 2)){
					continue;
				}
				delta = (1.0f/s->T[k] - 1.0f/s->T[k-1]) * (s->exE[s->rts[k-1]] - s->exE[s->rts[k]]);				
                double randme = pcg32randn(trng);
				//printf("\ndelta = %f e^-delta = %f (rand = %f)...", delta, exp((double)-delta), randme);
				//getchar();
				// do the exchange 
				//printf("pstep = %i         exchange %i --- %i\n", p, k, k-1);
				//getchar();
				//if( delta < 0.0 || pcg32randn(s->rng + tid) < exp(-delta) ){
				if( delta < 0.0 || randme < exp(-delta) ){
					/* temp swap function */
					//printf("<before> spins for %i and %i\n", k-1, k);
					//printspins(s->dlat[s->rts[k-1]], s->L, 20);
					//printspins(s->dlat[s->rts[k]], s->L, 20);
					/* swap temperatures */
                    //printf("  YES\n");
					swap(s, k-1, k);
					/* global counter */
					ex++;
					/* local counter */
					s->ex[k] += 1.0f;
					//printf("%i <--> %i\n", k, k-1);
					//printf("<after> spins for %i and %i\n", k-1, k);
					//printspins(s->dlat[s->m[k-1]], s->L, 20);
					//printspins(s->dlat[s->m[k]], s->L, 20);
					//printarray<int>(s->rts, s->R, "rts");
					//printarray<float>(s->T, s->R, "T");
					//printf("\n");
					//printarray<int>(array, s->R, "R");
					//printarray<int>(s->trs, s->R, "trs");
					//getchar();
				}
                else{
                    //printf("NO\n"); fflush(stdout);
                }
				/* set energy to zero */
				//s->exE[s->rts[k]] = 0.0;
		}
		//for(int i = 0; i < s->R; ++i){
		//	printf("exE[%i] = %f\n", i, s->exE[s->rts[i]]);
		//}	
		/* set energy of last replica to zero */
		//s->exE[s->rts[s->R-1]] = 0.0;
		//getchar();
	}
	/* sync again */
	#pragma omp barrier
	return ex;
}

/* measure phase occurs at mzone */
#ifdef MEASURE
void measure(setup_t *s, int tid, int a, int b, int p, pcg32_random_t *trng){
	if( p >= s->mzone ){
		for(int i=0; i< s->blocks; i++){
			reset_mcmc_statistics( s, tid, a, b);
			simulation(s, tid, a, b, trng); 
			accum_block_statistics( s, tid, a, b );
		}
	}
}
#endif

#endif
