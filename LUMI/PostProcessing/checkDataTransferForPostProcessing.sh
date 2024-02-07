#!/bin/bash -l
set -e

outputPath=

PDBids=($(ls ${outputPath}))

echo "Folders in ${outputPath} that miss the TAR archive with topology information:"
for PDBid in ${PDBids[@]}; do
  if ! [ -d ${outputPath}/${PDBid} ]
  then
    continue
  fi
  if ! [ -f ${outputPath}/${PDBid}/${PDBid}.tar.gz ]
  then
    echo ${PDBid}
  fi
done
