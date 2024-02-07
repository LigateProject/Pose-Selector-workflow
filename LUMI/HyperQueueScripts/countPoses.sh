#!/bin/bash
set -e

DATA=

success=0
failure=0

cd ${DATA}
for i in $(ls); do
  if (( $(ls ${i}/*.tar.gz | wc -w) > 0 ))
  then
    success=$(( ${success}+$(ls ${i} | grep -c "tar") ))
  fi
  if [ -f ${i}/failed.txt ]
  then
    failure=$(( ${failure}+$(cat ${i}/failed.txt | wc -l) ))
  fi 
done

echo "Successfully simulated poses:"
echo ${success}

echo "Failed poses:"
echo ${failure}
