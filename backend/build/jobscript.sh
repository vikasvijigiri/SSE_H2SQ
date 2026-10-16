#!/bin/bash
#PBS -o $PBS_O_WORKDIR/run.out
#PBS -e $PBS_O_WORKDIR/run.err
##PBS -q phi
#PBS -l nodes=1:ppn=1  # Request 1 node with 1 core

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > pbs

############################## COMPILATION & EXECUTION #####################################

##g++ -O3 -std=c++17 -Wall -pedantic-errors -g -msse2 -march=native -funroll-loops -ffast-math -fomit-frame-pointer -fstrict-aliasing -o ./main -ggdb ./*.cpp  
make
time ./sse_sim_test input_params.in  # Run in background, redirect output to log file
# End of jobscript.sh