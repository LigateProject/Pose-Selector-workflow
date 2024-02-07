#!/bin/bash

#SBATCH --job-name=PSFolderPreparation
#SBATCH --output=PSFolderPreparation.o%j
#SBATCH --error=PSFolderPreparation.e%j
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00:00

set -e

# So far, this script has not been used. The TPR files produced with prepareFolderStructure.sh work fine although they were generated with a GPU version of GROMACS 2023.2.

module use /appl/local/csc/modulefiles
module load gromacs/2023.2

PATH_TO_INPUT=

PDBids=($(ls ${PATH_TO_INPUT}))

for PDBid in ${PDBids[@]}; do
  # skip folders where we already started to simulate
  if (( $(find ${PATH_TO_INPUT}/${PDBid} -name *.tar.gz | wc -l)>0 ))
  then
    continue
  fi
  # skip if all poses in a folder could not be simulated
  if (( $(ls -I "fail*" ${PATH_TO_INPUT}/${PDBid} | wc -w)==0 ))
  then
    continue
  fi
  pose_list=($(ls -I "fail*" ${PATH_TO_INPUT}/${PDBid}))
  # we still need to set up the folder structure for a system if the last pose doesn't have folders for all replicas
  if ! [ -d ${PATH_TO_INPUT}/${PDBid}/${pose_list[-1]}/rep_8 ]
  then
    echo "Generating additional .tpr files for ${PDBid} ..."
    for pose in ${pose_list[@]}; do
      for i in {1..8}; do
        mkdir -p "$PATH_TO_INPUT/$PDBid/$pose/rep_$i"
      done
      mv "$PATH_TO_INPUT/$PDBid/$pose/simulation.tpr" "$PATH_TO_INPUT/$PDBid/$pose/rep_1/"
      for i in {2..8}; do
        gmx convert-tpr -generate_velocities -velocity_temp 298 -velocity_seed -1 -s "$PATH_TO_INPUT/$PDBid/$pose/rep_1/simulation.tpr" -o "$PATH_TO_INPUT/$PDBid/$pose/rep_$i/simulation.tpr" > "$PATH_TO_INPUT/$PDBid/$pose/rep_$i/seed_info.log"
      done	
    done
    echo "Done with ${PDBid}!"
  fi
done
