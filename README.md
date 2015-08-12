# mctrueke
multi-core version of trueke, a Monte Carlo simulator for the 3D random field Ising model

# compile
# note: it may be necessary to adapt the Makefile for your machine
make clean
make

# how to run
./bin/mctrueke -l <L> <R> -t <T> <dT> -h <h> -s <pts> <mzone> <drop> <mcs> <meas> <period>
-r <r> -x <rthreads> <sthreads>

# example
bin/mctrueke 128 3.0 0.1 0.0 16 10 2 2 2 2 1
