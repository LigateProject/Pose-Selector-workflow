#!/bin/bash
set -e

module use /appl/local/csc/modulefiles # for access to pre-installed GROMACS versions

module load gromacs/2023.2 # for consistency

module load cray-python/3.9.12.1

# additional Python packages
module load typer/0.9.0-cpu

# HQ itself
module load hq/0.17.0-cpu

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
