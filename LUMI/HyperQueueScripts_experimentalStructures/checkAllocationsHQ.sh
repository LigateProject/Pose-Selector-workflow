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

hq --server-dir=$(pwd)/hq-servers alloc list

# calculate number of allocations to be checked
hq --server-dir=$(pwd)/hq-servers alloc list >> HQStatus.txt
iterations=$(( $(cat HQStatus.txt | wc -l)-4 ))
rm HQStatus.txt

# display allocations one by one
count=4
for h in $(seq 1 ${iterations}); do
hq --server-dir=$(pwd)/hq-servers alloc list >> HQStatus.txt
number=$(tail -n +${count} HQStatus.txt | head -1 | awk '{print $2}')
rm HQStatus.txt
echo ${number}
hq --server-dir=$(pwd)/hq-servers alloc info ${number}
count=$(( ${count}+1 ))
done
