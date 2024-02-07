#!/bin/bash -l
set -e

notReady=()
ready=()
completed=()

for PDBid in $(ls); do
  if [ -f ${PDBid}/${PDBid}.tar.gz ]
  then
    if (( $(ls ${PDBid}/*.tar.gz | wc -w)==$(ls ${PDBid}/*.tar.gz | wc -w)+1 ))
    then
      completed+=(${PDBid})
    else
      ready+=(${PDBid})
    fi
  else
    notReady+=(${PDBid})
  fi
done

echo "Folders ready for post-processing:"
echo ${ready[@]}
echo "Folders ready for implicit-solvent calculations:"
echo ${completed[@]}
echo "Folders not yet ready for post-processing:"
echo ${notReady[@]}
