#!/bin/bash -l
set -e

module use /appl/local/csc/modulefiles
module load gromacs/2023.2

# set environment variables
NTMPI=16
NTOMP=1

export OMP_NUM_THREADS=${NTOMP}

if [[ -z "${NUM_ITERS}" ]]; then
  echo "NUM_ITERS has to be specified"
  exit 1
fi

# CPU-GPU mapping on LUMI
#    CPU cores      NUMA     socket
#     0-7,    8-15    0        0
#   16-23,   24-31    1        0
#   32-39,   40-47    2        0
#   48-55,   56-63    3        0
#   64-71,   72-79    0        1
#   80-87,   88-95    1        1
#  96-103, 104-111    2        1
# 112-119, 120-127    3        1

# different core numberings
# GROMACS pinning ID:   0    1    2    3
# OS core ID:           0  128    1  129
# hwloc:                0    0    1    1

# hwloc considers only physical cores, does not see cores blocked by the operating system and uses a consecutive numbering (not GROMACS's numbering)
HWLOCBEGIN=(dummy 0 16 32 48 64 80 96 112)
HWLOCEND=(dummy 15 31 47 63 79 95 111 127)
# idea: one simulation per NUMA node
PINOFFSET=(dummy 0 32 64 96 128 160 192 224)

# run simulations
pids=()

# double-check HyperQueue
echo ${CWD} $(date) $(hostname) >> HyperQueueRuns/currentlyRunning.txt

# work in the nodes memory /tmp
cd /tmp
cp -r ${CWD} .
archive=$(echo ${CWD} | rev | cut -d "/" -f 1 | rev)
cd ${archive}
# signal to analysis script that the job was started
mkdir -p ${CWD}/rep_1/cpt

for (( i=1; i<=${NUM_ITERS}; i++ )); do
    cd rep_${i}
    mkdir -p cpt

    OMP_PLACES=cores OMP_PROC_BIND=close hwloc-bind --cpubind "core:${HWLOCBEGIN[${i}]}-${HWLOCEND[${i}]}" \
    gmx mdrun  -pin on -pinoffset ${PINOFFSET[${i}]} -pinstride 2 \
               -deffnm simulation \
               -cpnum -cpt 60 -cpo cpt/state \
               -ntmpi ${NTMPI} -ntomp ${NTOMP} \
               -maxh 4 &

    # Store the PID for later
    pids+=($!)
    cd ..
done

# Wait for all pids
# This approach makes sure that wait will fail if any of the children have failed
for pid in ${pids[@]}; do
    wait $pid
done

# clean up to reduce the number of files (only for successfully run simulations)
if (( $(find rep_* -name simulation.gro | wc -l) == ${NUM_ITERS} ))
then
    # remove all GROMACS files but the .tpr and .xtc
    rm -r rep_*/cpt/ rep_*/simulation.edr rep_*/simulation.gro rep_*/simulation.log
    rm -r hq-work
    # HQ needs the output file to be located on the file system
    cp ${CWD}/hq-work/job-*/stderr .
    # convert pose into tar archive
    cd ..
    tar czf ${archive}.tar.gz ${archive}
    cp -r ${archive}.tar.gz ${CWD}/..
    rm -r ${CWD}
    rm -r ${archive} ${archive}.tar.gz
fi
