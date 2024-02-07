#!/bin/bash

#PBS -q qgpu_eurohpc
#PBS -N ASMESIS-SV_Validation
#PBS -l select=1:ngpus=8:mpiprocs=1:ompthreads=4,walltime=48:00:00

set -e; shopt -s expand_aliases

module load GROMACS

# node info (CPU only)
## physical id: 0/1 => sockets: 2
## cpu cores: 64 => 64 cores per socket and 128 cores in total
## siblings: 64 => no hyper-threading enabled
## processor: 0 ... 127 => 128 threads => 1 thread per core => no hyper-threading

# run simulation (not tested)
## GROMACS was compiled with GPU support (use later in dd-22-10) and for thread-MPI, not MPI
## GROMACS tries to be smart and guess the number of thread-MPI ranks and OpenMP threads by itself
## doubly indicating these numbers to be on the safe side
## play around with number of thread-MPI ranks and OpenMP threads for optimal performance
## for dd-22-10: a node with 8 GPUs can probably only be saturated with 8 simulations => get partial nodes or load 8 simulations on one node?

PATH_TO_SIM=(run_01 run_02 run_03 run_04 run_05 run_06 run_07 run_08)

tMPIranks=4 # this is not MPI, thus mpiprocs=1 (do not use MPI)
OMPthreads=4 # keep in sync with ompthreads option of PBS
offset=()
PIDs=()

for i in $(seq 0 7); do

value=$(($i * $tMPIranks * $OMPthreads))
offset+=($value)

cd $inputFileLocation
echo $inputFileLocation
cd ${PATH_TO_SIM[${i}]}
mkdir cpt
# When is it useful to have separate PME ranks?
# update gpu fails with some weird error message
gmx mdrun -deffnm simulation -pin on -cpnum -cpt 60 -cpo cpt/state -ntmpi ${tMPIranks} -ntomp ${OMPthreads} -pinstride 1 -pinoffset ${offset[${i}]} -gpu_id ${i} -maxh 48 -nb gpu -pme gpu -pmefft gpu -bonded gpu -update cpu -npme 1 &
PIDs+=($!)
cd ..

done

wait ${PIDs[0]} ${PIDs[1]} ${PIDs[2]} ${PIDs[3]} ${PIDs[4]} ${PIDs[5]} ${PIDs[6]} ${PIDs[7]}

