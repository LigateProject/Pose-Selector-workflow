#!/bin/bash

# Data paths
## local cluster
ROOT_PATH=${ligate}/PSLargeScale
PATH_FILES_TO_TRANSFER=${ROOT_PATH}/inputFiles
PATH_TRANSFER_METADATA=${ROOT_PATH}/transferredSystems

## external cluster
PATH_ON_LUMI=

# copy files
BATCH_SIZE=100

systemsToTransfer=($(cat ${PATH_TRANSFER_METADATA}/toBeTransferredToLumi.txt))
transferredSystems=($(cat ${PATH_TRANSFER_METADATA}/transferredToLumi.txt))
count=0
index=0
#start_time=$(date +%s.%N)
while (( ${count} < ${BATCH_SIZE} )); do
  if ! $(echo ${transferredSystems[@]} | grep -q "${systemsToTransfer[${index}]}")
  then
    rsync -aPhv ${PATH_FILES_TO_TRANSFER}/${systemsToTransfer[${index}]} lumi:${PATH_ON_LUMI} --exclude *.itp --exclude *.ndx --exclude *.top --exclude *.gro --exclude incorrect_number_of_disulfide_bond
    wait $!
    echo "All poses of ${systemsToTransfer[${index}]} have been copied!"
    count=$(( ${count}+1 ))
    echo ${systemsToTransfer[${index}]} >> ${PATH_TRANSFER_METADATA}/transferredToLumi.txt
  fi
  index=$(( ${index}+1 ))
done
#end_time=$(date +%s.%N)
#echo "$end_time - $start_time" | bc > execution_time.txt
echo "The end"
