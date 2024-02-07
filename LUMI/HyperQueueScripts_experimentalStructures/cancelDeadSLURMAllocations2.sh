#!/bin/bash
set -e

slurmID=""
user=""

DATA=

# get current list of nodes corresponding to SLURM allocations that are attributed to HyperQueue
slurmList=($(squeue -u ${user} -t RUNNING | tail -n +2 | grep "hq-alloc" | awk '{print $1}'))
slurmNodeList=($(squeue -u ${user} -t RUNNING | tail -n +2 | grep "hq-alloc" | awk '{print $8}'))

lostAllocations=""

for i in $(seq 0 $(( $(echo ${slurmList[@]} | wc -w)-1 ))); do
  occurrences=$(tail -$(( $(echo ${slurmList[@]} | wc -w) * 5 )) ${DATA}/HyperQueueRuns_experimentalStructures/currentlyRunning.txt | grep "${slurmNodeList[${i}]}") || lostAllocations="${lostAllocations} ${slurmList[${i}]}"
done

echo "Running allocations that did not contribute to the last $(( $(echo ${slurmList[@]} | wc -w) * 5 )) simulations (5 times the number of currently running SLURM allocations):"
echo ${lostAllocations}

for i in ${lostAllocations}; do
  scancel ${i}
done

echo "Cancelled all SLURM allocations listed in the line above."
