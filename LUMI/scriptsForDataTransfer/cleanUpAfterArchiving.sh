#!/bin/bash
set -e

ROOT_PATH=${ligate}/PSLargeScale

PDBids=$(ls ${ROOT_PATH}/archivedSystems | cut -d "." -f 1)

for PDBid in ${PDBids[@]}; do
  if [ -d ${ROOT_PATH}/inputFiles/${PDBid} ]
  then
    rm -r ${ROOT_PATH}/inputFiles/${PDBid}
  fi
done
