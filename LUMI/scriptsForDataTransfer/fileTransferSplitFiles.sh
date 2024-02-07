#!/bin/bash
set -e

# Data paths
## local cluster
ROOT_PATH=${ligate}/PSLargeScale
PATH_FILES_TO_TRANSFER=${ROOT_PATH}/inputFiles
PATH_TRANSFER_METADATA=${ROOT_PATH}/transferredSystems

# fraction of files that goes to LUMI,  others go to MeluXina
fractionToSendToLumi=0.7

# keep track of systems to be transferred
totalNumberOfSystems=$(ls ${PATH_FILES_TO_TRANSFER} | wc -w)
numberOfSystemsOnLumi=0 # $(cat ${PATH_TRANSFER_METADATA}/toBeTransferredToLumi.txt | wc -l) only has to be evaluated if systems have not yet been archived
numberOfSystemsOnMeluXina=$(cat ${PATH_TRANSFER_METADATA}/toBeTransferredToMeluXina.txt | wc -l)
numberToTransfer=$(( ${totalNumberOfSystems} - ${numberOfSystemsOnLumi} - ${numberOfSystemsOnMeluXina} ))

SYSTEMS_TO_TRANSFER=$(ls -ltr ${PATH_FILES_TO_TRANSFER} | tail -n +2 | rev | awk '{print $1}' | rev | tail -n ${numberToTransfer})

# create lists of PDB IDs to be copied to LUMI and MeluXina
numberToLumi=$(printf "%.0f" $(echo ${numberToTransfer} ${fractionToSendToLumi} | awk '{print $1 * $2}'))
## just copy the first ${numberToLumi} systems to LUMI, the rest is copied to MeluXina
count=0
for PDBid in ${SYSTEMS_TO_TRANSFER[@]}; do
  if (( ${count} < ${numberToLumi} ))
  then
    echo ${PDBid} >> ${PATH_TRANSFER_METADATA}/toBeTransferredToLumi.txt
  else
    echo ${PDBid} >> ${PATH_TRANSFER_METADATA}/toBeTransferredToMeluXina.txt
  fi
  count=$(( ${count} + 1 ))
done
