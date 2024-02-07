#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

PATH_TO_MDP=scriptsPS/mdp

PDBid=$(pwd | rev | cut -d "/" -f 1 | rev)

for pose in pose_*; do
# need to add ions to get a neutral simulation box
## It should be -maxwarn 1, but right now, we have some issues with STaGE not re-naming atoms properly (ugly fix that must be replaced)
gmx grompp -f ${PATH_TO_MDP}/em.mdp -c ${pose}/solvated.gro -r ${pose}/solvated.gro -p topol_amber_${pose}.top -o ${pose}/addIons.tpr -maxwarn 2 || true

## error catching
## addIons.tpr must exist
if ! [ -f ${pose}/addIons.tpr ]
then
echo "For ${pose} of ${PDBid}, ions could not be added to the complex structure (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

set +e
gmx genion -s ${pose}/addIons.tpr -o ${pose}/ions.gro -p topol_amber_${pose}.top -pname NA -nname CL -neutral <<-eof
SOL
eof
set -e

## error catching
## ions.gro must exist and must not be empty
if ! [ -s ${pose}/ions.gro ]
then
echo "For ${pose} of ${PDBid}, ions could not be added to the complex structure (genion error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

rm mdout.mdp ${pose}/addIons.tpr ${pose}/solvated.gro 
if [ -f "#topol_amber_${pose}.top.1#" ]
then
rm "#topol_amber_${pose}.top.1#"
echo "For ${pose} of ${PDBid}, ions were added. Removing topology back-up."
else
echo "For ${pose} of ${PDBid}, no ions were added."
fi

echo "For ${pose} of ${PDBid}, the workflow step 'adding ions' has been completed successfully!"

done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the workflow step 'adding ions' could not be completed successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
