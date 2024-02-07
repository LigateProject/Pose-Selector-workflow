#!/bin/bash
set -e

PDBids=($(ls PSWorkflow))

# make a rough estimate of the resources needed
export LC_NUMERIC="en_US.UTF-8"

timeLimit=12

echo "Estimate time demand to run all jobs scheduled:"
timeDemand=0
for PDBid in ${PDBids[@]}; do
  DATA=${PDBid}
  # skip PDBid if there are no more files or directories containing the pattern "pose"
  var=$(ls ${DATA}/pose*) || continue
  # only run simulations if there are files with the pattern "pose" in their name that have not been converted to archives yet
  if (( $(ls ${DATA}/pose* | wc -w)>$(ls ${DATA}/pose_*.tar.gz | wc -w) ))
  then
    numberOfPoses=$(ls ${DATA} | grep "pose" | grep -v "tar.gz" | wc -l)
    # one pose usually requires 2.5 min
    timeDemandPose=$(echo ${numberOfPoses} 24 | awk '{print $1/$2}')
    timeDemand=$(echo ${timeDemand} ${timeDemandPose} | awk '{print $1+$2}')
    echo "${PDBid} has ${numberOfPoses} poses => ${timeDemandPose} h"
  fi
done

echo "Total time demand: ${timeDemand} h"
if (( $(printf "%.0f\n" "${timeDemand}") == 0 ))
then
  echo "No jobs to be submitted! Stopping script."
  exit
fi
numberOfAllocations=$(echo ${timeDemand} ${timeLimit} | awk '{print $1/$2}')
# round to full integer
numberOfAllocations=$(echo ${numberOfAllocations} 1 | awk '{print $1+$2}')
numberOfAllocations=$(printf "%.0f\n" "${numberOfAllocations}")
echo "Requesting ${numberOfAllocations} nodes via HyperQueue."

# start HQ server
mkdir ../HyperQueueRuns/startHQ
cd ../HyperQueueRuns/startHQ
nohup bash ../../HyperQueueScripts/startHQ.sh >> nohup.out &

echo "Started HQ server"
sleep 10
grep "Pid" nohup.out

cp ../../HyperQueueScripts/createAllocationsHQ.sh .
sed -i createAllocationsHQ.sh -e "s/numberOfAllocations=0/numberOfAllocations="${numberOfAllocations}"/g"
sed -i createAllocationsHQ.sh -e "s/timeLimit=x/timeLimit="${timeLimit}"/g"
bash createAllocationsHQ.sh

cd ..

sleep 5

# do actual work
for PDBid in ${PDBids[@]}; do
  DATA=${PDBid}
  # skip PDBid if there are no more files or directories containing the pattern "pose"
  var=$(ls ${DATA}/pose*) || continue
  # only run simulations if there are files with the pattern "pose" in their name that have not been converted to archives yet
  if (( $(ls ${DATA}/pose* | wc -w)>$(ls ${DATA}/pose_*.tar.gz | wc -w) ))
  then
    echo "Submitting ${PDBid} ..."
    mkdir -p ${PDBid}
    cd ${PDBid}
    cp ../../HyperQueueScripts/runHQ.sh runHQ_${PDBid}.sh
    sed -i runHQ_${PDBid}.sh -e "s/PDBid=/PDBid=${PDBid}/g"
    cp ../../HyperQueueScripts/hq-md-run.py .
    nohup bash runHQ_${PDBid}.sh >> nohup.out &
    echo $! >> pid.txt
    cd ..
  fi
done

sleep 15

# start script that removes SLURM allocations that are no longer used by HyperQueue
nohup bash HyperQueueScripts/cleanUpSLURM.sh >> nohup.out &
echo $! >> pid.txt

# Don't forget to shut HQ server down when simulations are done
