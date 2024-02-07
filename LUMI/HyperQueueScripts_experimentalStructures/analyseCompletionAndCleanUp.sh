#!/bin/bash
set -e

DATA=

success=()
failure=()
pending=()

cd ${DATA}
for i in $(ls); do
  cd ${i}
  # If all existing directories have been converted into tar archives, we're fine
  if (( $(ls -I "fail*" | wc -w) == $(ls | grep -c "tar") ))
  then
    success+=(${i})
  else
    failed=0
    notSimulated=0
    for j in $(ls -I "*.tar.gz"); do
      # If a cpt folder is found, the simulation was submitted but couldn't be converted into a tar file due to some error with the MD simulation.
      if [ -d ${j}/rep_1/cpt ]
      then
        if ! [ -f failed.txt ]
        then
          touch failed.txt
        fi
        mkdir -p failed
        failed=$(( ${failed}+1 ))
	echo ${j} >> failed.txt
	mv ${j}/hq-work/job-*/stderr failed/stderr_${j}
	rm -r ${j}
      # The MD script first creates a cpt folder. If this folder can't be found, the simulation batch was not submitted.
      else
	notSimulated=$(( ${notSimulated}+1 ))
      fi
    done
    if (( ${failed} > 0 ))
    then
      failure+=(${i})
    fi
    if (( ${notSimulated} > 0 ))
    then
      pending+=(${i})
    fi
  fi
  cd ..
done

echo "MD simulations are complete for these complexes:"
echo ${success[@]}

echo "MD simulations failed for poses of these complexes:"
echo ${failure[@]}

echo "MD simulations have not been started for poses of these complexes:"
echo ${pending[@]}
