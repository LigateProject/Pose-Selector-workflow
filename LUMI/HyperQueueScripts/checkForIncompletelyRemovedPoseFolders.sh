#!/bin/bash -l
set -e

outputPath=

PDBids=($(ls ${outputPath}))

echo "Poses that were simulated but not removed:"
for PDBid in ${PDBids[@]}; do
  poses=$(ls ${outputPath}/${PDBid} -I "*.tar.gz" -I "fail*")
  for pose in ${poses}; do
    if [ -f ${outputPath}/${PDBid}/${pose}.tar.gz ]
    then
      echo ${PDBid}/${pose}
      ls ${outputPath}/${PDBid}/${pose}
    fi
  done
done
