#!/bin/bash

#SBATCH --job-name=PSDataArchiving
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH -t 48:00:00
#SBATCH -o slurm-%J.out

set -e

# Data paths
ROOT_PATH=${ligate}/PSLargeScale/experimentalStructures
PATH_FILES_TO_ARCHIVE=${ROOT_PATH}/inputFiles
PATH_OUTPUT=${ROOT_PATH}/archivedSystems
PATH_TRANSFER_METADATA=${ROOT_PATH}/transferredSystems

PDBids=($(cat ${PATH_TRANSFER_METADATA}/transferredToLumi.txt)) # $(cat ${PATH_TRANSFER_METADATA}/transferredToMeluXina.txt))

cwd=/scratch

mkdir -p ${cwd}
cd ${cwd}
for PDBid in ${PDBids[@]}; do
  # check if archive already exists
  if ! [ -e "${PATH_OUTPUT}/$PDBid.tar.gz" ]; then
    cp -r  ${PATH_FILES_TO_ARCHIVE}/${PDBid} .
    tar -czf ${PDBid}.tar.gz ${PDBid}
    cp ${PDBid}.tar.gz ${PATH_OUTPUT}
    echo "Archived ${PDBid}!"
    rm -r ${PDBid}.tar.gz ${PDBid}
  fi
done
echo "Done!"
