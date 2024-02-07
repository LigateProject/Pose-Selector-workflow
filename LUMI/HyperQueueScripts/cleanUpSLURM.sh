#!/bin/bash
set -e

# stay active in the background for 24 hours
for i in $(seq 1 96); do
  # stop if there are no more HyperQueue jobs in the queue
  reasonsToStayActive=$(squeue -u | grep "hq-alloc") || exit
  echo "Starting control cycle ${i}"
  # sleep the quarter of an hour
  sleep 900
  # clean up SLURM allocations
  bash HyperQueueScripts/cancelDeadSLURMAllocations2.sh
done
