#!/bin/bash
set -e; shopt -s expand_aliases

PDBids=(1ajp 1ps3 1uto 1x8d 2r58 3dyo 3gv9 3i3b 4ahs 4jsz 4mre 4q90 4qsu 4u54 4y0a 5aol 5mz8 5z5f 6abx 6r9u)

for PDBid in ${PDBids[@]}; do
cd ${PDBid}
qsub -v inputFileLocation=$(pwd) scripts/runSimulationInQueueRestart.sh
done
