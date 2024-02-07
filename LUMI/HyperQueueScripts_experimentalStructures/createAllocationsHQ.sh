#!/bin/bash
set -e

# modules consistent with GROMACS
module load LUMI/22.08
module unload cray-mpich/8.1.27 # due to tMPI
module unload cray-libsci/22.08.1.1 # due to tMPI
module load buildtools/22.08 # loads CMake/3.24.0
module load cray-python/3.9.12.1

# additional Python packages
module load typer/0.9.0

# HQ itself
module load hq/0.16.0

# let HyperQueue create as many allocations as needed
numberOfAllocations=0
timeLimit=x

for i in $(seq 1 ${numberOfAllocations}); do
hq --server-dir=$(pwd)/hq-servers alloc add slurm --time-limit ${timeLimit}h --workers-per-alloc 1 --idle-timeout 60s --backlog 1 -- --account= --partition=standard-g
# additional arguments for small-g: --ntasks-per-node=8 --cpus-per-task=8 --gpus-per-node=8
done
