#!/bin/bash
set -e

BATCH_SIZE=101

sbatch structuresForMLModelTraining_experimentalStructures.sh
for i in $(seq 0 101 2323); do
  sed -i structuresForMLModelTraining_experimentalStructures.sh -e "s/START_INDEX="${i}"/START_INDEX="$(( ${i}+${BATCH_SIZE} ))"/g"
  sbatch structuresForMLModelTraining_experimentalStructures.sh
done
