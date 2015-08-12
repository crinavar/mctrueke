// -------------------------------------------------------------//
// pt: parallel tempering on the GPU                            //
// authors: Cristobal A. Navarro, Huang Wei, Youjin Deng.       //
// last modified: 6 Nov, 2014.                                  //
//--------------------------------------------------------------//
#ifndef _SETUP_H_
#define _SETUP_H_

/* function declarations */
void init(setup_t *s, int argc, int argv);
void printparams(setup_t *s);
void getparams(setup_t *s, int argc, char **argv);
void newseed(int *seed);
void malloc_arrays(setup_t *s);
void adjustparams(setup_t *s);

void adjustparams(setup_t *s){
    s->blocks = 1;
    /* total number of spins per replica */
    s->N = (s->L)*(s->L)*(s->L);
    /* adjust R to a multiple of threads; R' = ceil(R/threads) *threads */
    s->R = (int)ceil((float)s->R/(float)s->rthreads) * s->rthreads;
    /* measure zone */
    if( s->mzone == -1 ){
        s->mzone = (int) ((double)s->pts / log((double)s->pts / (double)s->L) );
    }
}

/* init */
void init(setup_t *s, int argc, char **argv){
	//printf("init....{\n"); fflush(stdout);
	/* get parameters */
	getparams(s, argc, argv);
    /* adjust some parameters related to memory pool and active replicas*/
    adjustparams(s);
#ifdef MEASURE
	/* folders for output */
	s->obsfolder = "data";
	s->plotfolder = "plots";
	make_output_folders(s->obsfolder, s->plotfolder);
#endif
	/* random seeds */
	s->seed = time(NULL);
	s->nseed = s->seed + 7919;
	/* fixed seeds */
	//s->seed = s->nseed = 7919;

	srand(s->nseed);
	/* set the number of threads as the number of GPUs */
	//omp_set_num_threads(s->rthreads);
    // allow nested parallelism
    omp_set_nested(1);
    // disable adjustment on the number of threads in nested parallel regions
    omp_set_dynamic(0);

	/* alocate main arrays */
	malloc_arrays(s);

	/* reset table of obersvables per realization */
	reset_realization_statistics(s, s->R);

	/* print parameters */
	printparams(s);
    //printarray<float>(s->T, s->R, "T");
    //printarray<float>(s->E, s->R, "E");
	//printf("}:ok\n\n"); fflush(stdout);
}

/* malloc arrays */
void malloc_arrays( setup_t *s ){
	/* allocate the main arrays */
	s->hlat 	= (int **)malloc(sizeof(int *) * s->R);

	/* T is a sorted temp array */
	s->T = (float*)malloc(sizeof(float)*s->R);
	/* ex is a per temperature counter array */
	s->ex = (float*)malloc(sizeof(float)*s->R);
	/* avex is a per temperature counter array */
	s->avex = (float*)malloc(sizeof(float)*s->R);
	/* index arrays */ 
	s->rts = (int*)malloc(sizeof(int)*s->R);
	s->trs = (int*)malloc(sizeof(int)*s->R);

	/* host values for each replica */
	s->E = (float*)malloc(sizeof(float)*s->R);
	s->exE = (float*)malloc(sizeof(float) * s->R);
	s->M = (int*)malloc(sizeof(int)*s->R);
	s->F1 = (float3*)malloc(sizeof(float3)*s->R);
	s->F2 = (float3*)malloc(sizeof(float3)*s->R);
	/* observables table */
	s->obstable = (obset_t*)malloc(sizeof(obset_t)*s->R);
	// memory for H array
	s->hH = (int*)malloc(sizeof(int) * s->N);
    // distributed rng 
    s->rng = (pcg32_random_t **)malloc(sizeof(pcg32_random_t*) * s->rthreads);

    
	/* multi-core setup -- replica level paralelism */
    const int prime = s->nseed + 17;
	#pragma omp parallel default(none) shared(s) proc_bind(spread) num_threads(s->rthreads)
	{
		int tid, nt, r, k;
		/* set threads */
		threadset(s, &tid, &nt, &r);
		/* malloc the data for 'r' replicas on each GPU */
		for(int j = 0; j < r; ++j){
			k = tid * r + j;
            /* replica allocation */
            s->hlat[k]= (int*)malloc(sizeof(int) * s->N);

            /* array of temperatures increasing order */
            s->T[k] = s->TR - (s->R-1 - k)*s->dT;

            /* exchange counters initialization */
            // warning: possible false sharing if R <= threads * 64bytes
            s->ex[k] = 0;

            /* initialize index arrays */
            s->rts[k] = s->trs[k] = k;
		}	
        s->rng[tid] = (pcg32_random_t*)malloc(sizeof(pcg32_random_t) * s->sthreads);
        #pragma omp parallel default(none) shared(s, tid) proc_bind(close) num_threads(s->sthreads) 
        {
            #pragma omp for
            for(int q = 0; q < s->sthreads; ++q){
                pcg32_srandom_r(&(s->rng[tid][q]), prime + tid*s->sthreads + q, prime + 19 + tid*s->sthreads + q);
            }
        }
	}
	/* host memory setup for each replica */
    printarray<float>(s->T, s->R, "T");
    printf("\n");
}

/* print parameters */
void printparams(setup_t *s){
	printf("\tparameters:{\n");
	printf("\t\tL:                            %i\n", s->L);
	printf("\t\tvolume:                       %i\n", s->N);
	printf("\t\t[TR,dT]:                      [%f, %f]\n", s->TR, s->dT);
	printf("\t\tmag_field h:                  %f\n", s->h);
	printf("\t\treplicas:                     %i\n", s->R);
	printf("\t\tptsteps:                      %i\n", s->pts);
	printf("\t\tmzone:                        %i\n", s->mzone);
	printf("\t\tdrop_steps:                   %i\n", s->ds);
	printf("\t\tmcsteps:                      %i\n", s->ms);
	printf("\t\tmeasure:                      %i\n", s->fs);
	printf("\t\tperiod:                       %i\n", s->period);
	printf("\t\tnblocks:                      %i\n", s->blocks);
	printf("\t\trealizations:                 %i\n", s->realizations);
	printf("\t\tseed:                         %i\n", s->seed);
	printf("\t\trthreads:                     %i\n", s->rthreads);
	printf("\t\tsthreads:                     %i\n\t}\n", s->sthreads);
}

/* get parameters */
void getparams(setup_t *s, int argc, char **argv){
	/* if the number or arguments is not correct, stop the program */
	if(argc != 21){
		printf("run as:\n./bin/mctrueke -l <L> <R> -t <T> <dT> -h <h> -s <pts> <mzone> <drop> <mcs> <meas> <period> -r <r> -x <rthreads> <sthreads>\n");
		exit(1);
	}
	else{
		for(int i=0; i<argc; i++){
			/* lattice size and number of replicas */
			if(strcmp(argv[i],"-l") == 0){
				s->L = atoi(argv[i+1]);	
				s->R = atoi(argv[i+2]);
			}
			/* get TR and dT */
			else if(strcmp(argv[i],"-t") == 0){
				s->TR = atof(argv[i+1]);
				s->dT = atof(argv[i+2]);
			}
			/* the magnetic field constant */
			else if(strcmp(argv[i],"-h") == 0){
				s->h = atof(argv[i+1]);
			}
			/* ptsteps, drop steps, mc steps, final steps */
			else if(strcmp(argv[i],"-s") == 0){
				s->pts = atof(argv[i+1]);
                s->mzone = atoi(argv[i+2]);
				s->ds = atof(argv[i+3]);
				s->ms = atof(argv[i+4]);
				s->fs = atof(argv[i+5]);
				s->period = atof(argv[i+6]);
			}
			/* number of measure blocks and realizations */
			else if(strcmp(argv[i],"-r") == 0){
				s->realizations = atof(argv[i+1]);
			}	
			/* number of gpus */
			else if(strcmp(argv[i],"-x") == 0){
				s->rthreads = atoi(argv[i+1]);
				s->sthreads = atoi(argv[i+2]);
			}
		}
	}
	if( (s->L % s->rthreads) != 0 ){
		fprintf(stderr, "lattice dimensional size must be multiples of %i", s->rthreads);
    }
}

/* generate new seed */
void newseed(setup_t* s){
	//s->nseed = time(NULL);
}
#endif
