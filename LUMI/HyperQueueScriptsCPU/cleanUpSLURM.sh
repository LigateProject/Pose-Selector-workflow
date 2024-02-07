#!/bin/bash
set -e

# stay active in the background for 1 week
for i in $(seq 1 672); do
  # stop if there are no more HyperQueue jobs in the queue
  reasonsToStayActive=$(squeue -u | grep "hq-alloc") || exit
  echo "Starting control cycle ${i}"
  # sleep the quarter of an hour
  sleep 900
  # clean up SLURM allocations
  bash HyperQueueScriptsCPU/cancelDeadSLURMAllocations2.sh
done
