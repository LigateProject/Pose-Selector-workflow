#!/bin/bash
set -e

BATCH_SIZE=101

sbatch finalArchiving_experimentalStructures.sh
for i in $(seq 0 101 2323); do
  sed -i finalArchiving_experimentalStructures.sh -e "s/START_INDEX="${i}"/START_INDEX="$(( ${i}+${BATCH_SIZE} ))"/g"
  sbatch finalArchiving_experimentalStructures.sh
done
