# mctrueke
multi-core version of trueke, a Monte Carlo simulator for the 3D random field Ising model

# compile
# note: it may be necessary to adapt the Makefile for your machine
make clean
make

# how to run
./bin/mctrueke -l <L> <R> -t <T> <dT> -h <h> -s <pts> <mzone> <drop> <mcs> <meas> <period> -r <r> -x <rthreads> <sthreads>

# example
bin/mctrueke -l 32 50 -t 4.7 0.02 -h 1.0 -s 2000 1000 10 1 1 -r 100 -x 2 1
