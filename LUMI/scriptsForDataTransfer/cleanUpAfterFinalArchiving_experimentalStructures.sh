#!/bin/bash
set -e

# assemble a list of tar archives that have been prepared for long-term storage
cd $pslt/experimentalStructures/simulationInputFiles
ARCHIVES=($(ls))

# remove corresponding tar archives in LIGATE_AWH_Benchmark
cd ${ligate}/PSLargeScale/experimentalStructures/archivedSystems

for archive in ${ARCHIVES[@]}; do
  if [ -f ${archive} ]
  then
    rm ${archive}
  fi
done
