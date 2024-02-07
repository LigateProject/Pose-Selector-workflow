#!/bin/bash
set -e

BATCH_SIZE=10

sbatch finalArchiving.sh
for i in $(seq 1 10 2581); do
  sed -i finalArchiving.sh -e "s/START_INDEX="${i}"/START_INDEX="$(( ${i}+${BATCH_SIZE} ))"/g"
  sbatch finalArchiving.sh
done
