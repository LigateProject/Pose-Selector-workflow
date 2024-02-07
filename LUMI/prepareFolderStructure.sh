#!/bin/bash

#SBATCH --job-name=PSFolderPreparation
#SBATCH --output=PSFolderPreparation.o%j
#SBATCH --error=PSFolderPreparation.e%j
#SBATCH --partition=small-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00:00

version=2023.2

set -e

# hwloc: /usr/bin/hwloc-*; version 2.6.0

module use /appl/local/csc/modulefiles # for semi-official hipSYCL

module load LUMI/22.08
module load hipsycl/0.9.4 # loads its own rocm
module unload cray-mpich/8.1.27 # due to tMPI
module unload cray-libsci/22.08.1.1 # due to tMPI
module load buildtools/22.08 # loads CMake/3.24.0
module load cray-fftw/3.3.10.1
module load cray-python/3.9.12.1

export LD_LIBRARY_PATH=/appl/lumi/SW/LUMI-22.08/G/EB/rocm/5.3.3/llvm/lib:$LD_LIBRARY_PATH # to be on the safe side: it's not in the module environment but needed by the clang compiler

module load gromacs/2023.2-sycl-tmpi

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
