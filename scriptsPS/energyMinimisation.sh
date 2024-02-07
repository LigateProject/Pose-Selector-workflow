#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

PATH_TO_MDP=scriptsPS/mdp

PDBid=$(pwd | rev | cut -d "/" -f 1 | rev)

for pose in pose_*; do
# run grompp to get input .tpr file
gmx grompp -f ${PATH_TO_MDP}/em.mdp -c ${pose}/ions.gro -r ${pose}/ions.gro -p topol_amber_${pose}.top -o ${pose}/EM.tpr -po ${pose}/EMout.mdp || true

## error catching
## EM.tpr must exist
if ! [ -f ${pose}/EM.tpr ]
then
echo "For ${pose} of ${PDBid}, the energy minimisation failed (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

# run energy minimisation
cd ${pose}
gmx mdrun -v -deffnm EM || true
cd ..

## error catching
## EM.gro must exist and must not be empty
if ! [ -s ${pose}/EM.gro ]
then
echo "For ${pose} of ${PDBid}, the energy minimisation failed (mdrun error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

# file clean up such that subsequent scripts can be used without any problems
mv ${pose}/EM.gro ${pose}/ions.gro
rm ${pose}/EM*

echo "For ${pose} of ${PDBid}, the energy minimisation has been completed successfully."

done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the energy minimisation could not be completed successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
