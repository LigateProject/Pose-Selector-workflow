#!/bin/bash

#SBATCH --job-name=PreparePSDataForLongTermStorage
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH -t 24:00:00
#SBATCH -o slurm-%J.out

set -e

# Paths
cwd=/scratch
PATH_ARCHIVES=${ligate}/PSLargeScale/experimentalStructures/archivedSystems

# Data
cd ${PATH_ARCHIVES}
ARCHIVES=($(ls))

# Select batch of archives to work on
BATCH_SIZE=101
START_INDEX=0 # start index for this batch

# Go to the node's local hard drive not to overload the file system
mkdir -p ${cwd}
cd ${cwd}

# Work on archives
index=${START_INDEX}

while (( ${index} < $(( ${BATCH_SIZE}+${START_INDEX} )) )); do

  if (( ${index}==${#ARCHIVES[@]} ))
  then
    break
  fi

  echo "Starting adding velocity seeds to ${ARCHIVES[${index}]}!"

  # copy relevant data
  cp ${PATH_ARCHIVES}/${ARCHIVES[${index}]} .
  folder=$(echo ${ARCHIVES[${index}]} | cut -d "." -f 1)
  if [ -d ${ps1}/experimentalStructures/${folder} ]
  then
    directory=${ps1}/experimentalStructures
  elif [ -d ${ps2}/experimentalStructures/${folder} ]
  then
    directory=${ps2}/experimentalStructures
  elif [ -d ${ps3}/experimentalStructures/${folder} ]
  then
    directory=${ps3}/experimentalStructures
  fi

  # inflate
  tar xzf ${ARCHIVES[${index}]}
  rm ${ARCHIVES[${index}]}

  # add log files with velocity seeds from raw trajectories to existing archives
  mkdir ${folder}_raw_trajectories
  cd ${folder}_raw_trajectories
  for pose in $(ls ${directory}/${folder} | grep "pose"); do
    cp ${directory}/${folder}/${pose} .
    tar xzf ${pose}
    pose_folder=$(echo ${pose} | cut -d "." -f 1)
    for i in $(seq 2 8); do
      cp ${pose_folder}/rep_${i}/seed_info.log ../${folder}/${pose_folder}/seed_info_rep_${i}.log
    done
    rm -r ${pose_folder} ${pose}
  done
  cd ..

  # compress
  tar czf ${ARCHIVES[${index}]} ${folder}

  # copy archive to final destination
  cp ${ARCHIVES[${index}]} $pslt/experimentalStructures/simulationInputFiles

  # clean up
  rm -r ${folder}_raw_trajectories
  rm -r ${folder}
  rm ${ARCHIVES[${index}]}

  echo "Done adding velocity seeds to ${ARCHIVES[${index}]}!"
  index=$(( ${index}+1 ))

done

# Completed
echo "Done!"
