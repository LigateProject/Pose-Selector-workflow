#!/bin/bash
set -e

BATCH_SIZE=10

sbatch structuresForMLModelTraining.sh
for i in $(seq 1 10 2581); do
  sed -i structuresForMLModelTraining.sh -e "s/START_INDEX="${i}"/START_INDEX="$(( ${i}+${BATCH_SIZE} ))"/g"
  sbatch structuresForMLModelTraining.sh
done
