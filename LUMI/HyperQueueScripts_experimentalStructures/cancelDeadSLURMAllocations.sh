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

# calculate number of HyperQueue allocations to be checked
hq --server-dir=$(pwd)/hq-servers alloc list >> HQStatus.txt
iterations=$(( $(cat HQStatus.txt | wc -l)-4 ))
rm HQStatus.txt

slurmID=""
user=""

for h in $(seq 1 ${iterations}); do
  # identify oldest allocation
  hq --server-dir=$(pwd)/hq-servers alloc list >> HQStatus.txt
  number=$(tail -n +4 HQStatus.txt | head -1 | awk '{print $2}')
  rm HQStatus.txt

  if ! [ -z ${number} ]
  then
    # get slurm id
    hq --server-dir=$(pwd)/hq-servers alloc info ${number} >> HQStatus.txt
    slurmID=$(echo ${slurmID})" "$(grep "RUNNING" HQStatus.txt | awk '{print $2}')
    slurmID=$(echo ${slurmID})" "$(grep "FAILED" HQStatus.txt | awk '{print $2}')
    rm HQStatus.txt
  fi
done

# get current list of SLURM allocations that are attributed to HyperQueue
slurmList=$(squeue -u ${user} -t RUNNING | tail -n +2 | grep "hq-alloc" | awk '{print $1}')
echo "List of SLURM allocations that HyperQueue seems to have forgotten about:"
echo ${slurmList}
# remove SLURM allocations that HyperQueue forgot about
for i in ${slurmList}; do
  echo ${slurmID} | grep ${i} || scancel ${i}
done
