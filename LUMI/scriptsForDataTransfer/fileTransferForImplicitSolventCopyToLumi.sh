#!/bin/bash
set -e

# Data paths
## local cluster
ROOT_PATH=${ligate}/PSLargeScale
PATH_FILES_TO_TRANSFER=${ROOT_PATH}/inputFiles
PATH_TRANSFER_METADATA=${ROOT_PATH}/transferredSystems

## external cluster
PATH_ON_LUMI=

# copy files
START_INDEX=2601
BATCH_SIZE=

systemsToTransfer=($(cat ${PATH_TRANSFER_METADATA}/toBeTransferredToLumi.txt))
for index in $(seq ${START_INDEX} $(( ${START_INDEX}+${BATCH_SIZE}-1 ))); do
  system=${systemsToTransfer[${index}]}
  mkdir ${system}
  cd ${system}
  # copy topologies
  pose=$(ls ${PATH_FILES_TO_TRANSFER}/${system}/topol_amber_pose_* | head -1 | rev | cut -d "_" -f 1 | rev | cut -d "." -f 1)
  cp ${PATH_FILES_TO_TRANSFER}/${system}/*.itp ${PATH_FILES_TO_TRANSFER}/${system}/topol_amber_pose_${pose}.top .
  # remove pose-specific info about number of water molecules and Na/Cl ions in solution
  head -$(( $(grep -n "SOL" topol_amber_pose_${pose}.top | tail -1 | cut -d ":" -f 1)-1 )) topol_amber_pose_${pose}.top >> topol_amber.top
  rm topol_amber_pose_${pose}.top
  # copy .gro files with starting structure
  poses=($(ls ${PATH_FILES_TO_TRANSFER}/${system}/topol_amber_pose_*))
  for pose in ${poses[@]}; do
    number=$(echo ${pose} | rev | cut -d "_" -f 1 | rev | cut -d "." -f 1)
    # remove water molecules from .gro file
    head -$(( $(grep -n "SOL" ${PATH_FILES_TO_TRANSFER}/${system}/pose_${number}/ions.gro | head -1 | cut -d ":" -f 1)-1 )) ${PATH_FILES_TO_TRANSFER}/${system}/pose_${number}/ions.gro >> ions_pose_${number}.gro
    oldAtomCount=$(head -2 ions_pose_${number}.gro | tail -1)
    newAtomCount=$(( $(cat ions_pose_${number}.gro | wc -l)-2 ))
    sed -i ions_pose_${number}.gro -e "s/${oldAtomCount}/${newAtomCount}/g"
    tail -1 ${PATH_FILES_TO_TRANSFER}/${system}/pose_${number}/ions.gro >> ions_pose_${number}.gro
  done
  cd ..
  # store in archive, transfer to LUMI and clean up
  tar czf ${system}.tar.gz ${system}
  rsync -aPhv ${system}.tar.gz lumi:${PATH_ON_LUMI}/${system}
  rm -r ${system}.tar.gz ${system}
done
