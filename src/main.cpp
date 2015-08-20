// -------------------------------------------------------------//
// pt: parallel tempering on the GPU                            //
// authors: Cristobal A. Navarro, Huang Wei, Youjin Deng.       //
// last modified: 6 Nov, 2014.                                  //
//--------------------------------------------------------------//

/* comment MEASURE if you just want to check simulation performance. */
#define MEASURE

/* warning: changing the block dimensions can lead to undefined behavior. */
#define WARPSIZE 32
#define BX	32
#define BY	8
#define	BZ	4
#define BLOCKSIZE1D 512
#define Q_DIST 2000
#define PI 3.14159265359

// openmp defines
#define CHUNCKSIZE 128

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <cstring>
#include <omp.h>

#include <sys/time.h>
#include <syslog.h>

// gperftools
#include <gperftools/profiler.h>

/* local includes */
#include "pcg_basic.h"
#include "structs.h"
#include "cputools.h"
#include "tools.h"
#include "setup.h"
#include "pt.h"

int main(int argc, char **argv){
    //printf("uint64_t = %i bytes              int = %i bytes\n", sizeof(uint64_t), sizeof(int));
	/* openmp thread values */
	int tid, nt, r, a, b; 
	printf("\n**************** mctrueke ****************\n\n");
	/* setup handles the variables */
	setup_t s;

	/* initialization takes care of memory allocation */
	init(&s, argc, argv);	
	/* measure time */
    double t1, t2;
    t1 = omp_get_wtime();
	/* main simulation */ 
	for(int i = 0; i < s.realizations; i++){
		printf("[realization %i of %i]\n", i+1, s.realizations); fflush(stdout);
		/* multi-core PT simulation */
        pcg32_random_t trng;
		#pragma omp parallel default(none) shared(s) private(tid, nt, r, a, b, trng) proc_bind(spread) num_threads(s.rthreads)
		{
            double st1, st2;
			/* set the thread */
			threadset(&s, &tid, &nt, &r);
			a = tid * r;
			b = a + r;
			/* reset some data at each realization*/
			reset(&s, tid, a, b);
			
			/* distribution for H */
			hdist(s.hH, s.N, tid, a, b, &trng);

			/* equilibration */
            if(tid == 0){
                st1 = omp_get_wtime();
            }
            // equilibration
			metropolis(&s, tid, a, b, s.ds, &trng);
            #pragma omp barrier
            if(tid == 0){
                st2 = omp_get_wtime();
                printf("spin time = %f secs\n", st2-st1);
            }
			/* parallel tempering */
			pt(&s, tid, a, b, &trng);
			#ifdef MEASURE
				/* accumulate realization statistics */
				accum_realization_statistics( &s, tid, a, b, s.realizations );
			#endif
		}
		/* try a new seed */
		newseed(&s);
	}
    t2 = omp_get_wtime();

#ifdef MEASURE
	/* write physical results */
	physical_results(&s);
#endif
	/* total time */
	printf("ok: total time %.2f secs\n", t2-t1);
	printf("spinperf: %f ns\n", (t2-t1)/((double)s.R * (double)s.N) * (1000000000.0 / (double)s.ds));

	/* free memory */
	freemem(&s);
}	
