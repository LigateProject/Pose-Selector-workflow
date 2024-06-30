#!/bin/bash

#PBS -q qcpu_eurohpc
#PBS -N ASMESIS-SV_Validation
#PBS -l select=1:mpiprocs=1:ompthreads=1,walltime=48:00:00

set -e; shopt -s expand_aliases

module load GROMACS
module load Amber
module load parmEd/3.4
module load gmx_MMPBSA/1.5.6

source GMXMMPBSA.sh

PATH_TO_SCRIPTS=scriptsPS
LIGANDFF="gaff2"

solvationModel="gb" # Generalised Born model of implicit solvation
# a full list of input parameters for 3D-RISM and descriptions can be found at https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/#rism-namelist-variables
# a full list of input parameters for Poisson-Boltzmann and descriptions can be found at https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/#pb-namelist-variables
# a full list of requirements to use gmx_MMPBSA for single protein-ligand trajectories: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Protein_ligand/ST/

#PDBids=(6mj7 5yj8 5ovp 4mn3 5wxh 6j3p 4wn5 4zyf 6udu 6udt 5d3c 6dif 5dgw 6e9a)
PDBids=(6mj7)

for PDBid in ${PDBids[@]}; do
for i in $(seq -w 01 08); do
cd ${PDBid}/run_${i}

# grompp (unfortunately, the GROMACS versions on our cluster and on Karolina aren't the same such that TPR files have to be re-generated)
gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/simulationAnalysis.mdp -c rerunTPR.gro -p ../topol_rerun.top -n indexRerun.ndx -o rerun.tpr -po rerunOut.mdp -maxwarn 1

# filter trajectory
gmx trjconv -f ${solvationModel}.xtc -s rerun.tpr -n indexRerun.ndx -o ${solvationModel}_filtered.xtc -skip 20 <<-eof
complex
complex
eof

# run gmx_MMPBSA (see gmx_MMPBSA -h for meaning of flags; cannot open GUI part of gmx_MMPBSA => suppress corresponding error)
gmx_MMPBSA -O -i ${PATH_TO_SCRIPTS}/mmpbsa_99sbildn_${LIGANDFF}_${solvationModel}.in -cs rerun.tpr -cp ../topol_rerun.top -ci indexRerun.ndx -cg 1 0 -ct ${solvationModel}_filtered.xtc -o ${solvationModel}_mmpbsa.dat -eo ${solvationModel}_mmpbsa.csv || true &
PIDs+=($!)
sleep 60

done
done

wait ${PIDs[0]} ${PIDs[1]} ${PIDs[2]} ${PIDs[3]} ${PIDs[4]} ${PIDs[5]} ${PIDs[6]} ${PIDs[7]}
