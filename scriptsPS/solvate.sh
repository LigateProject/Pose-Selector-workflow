#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

PDBid=$(pwd | rev | cut -d "/" -f 1 | rev)

for pose in pose_*; do
# place the complex in a rhombic dodecahedron of the desired size (1.5 nm distance to the edges)
gmx editconf -f ${pose}/full.gro -o ${pose}/correctBox.gro -bt dodecahedron -d 1.5 || true

## error catching
## correctBox.gro must exist and must not be empty
if ! [ -s ${pose}/correctBox.gro ]
then
echo "For ${pose} of ${PDBid}, the complex structure could not be placed in a box of the correct size. Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose}
continue
fi

# solvate with TIP3P water
cp topol_amber.top topol_amber_${pose}.top
gmx solvate -cp ${pose}/correctBox.gro -cs spc216.gro -p topol_amber_${pose}.top -o ${pose}/solvated.gro || true

## error catching
## solvated.gro must exist and must not be empty
if ! [ -s ${pose}/solvated.gro ]
then
echo "For ${pose} of ${PDBid}, the complex structure could not be solvated (.gro file not obtained). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
## water molecules must be added to the topology
elif ! $(diff 'topol_amber_'${pose}'.top' '#topol_amber_'${pose}'.top.1#' | grep -q "SOL")
then
echo "For ${pose} of ${PDBid}, the complex structure could not be solvated (water molecules were not added to topology). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

rm "#topol_amber_${pose}.top.1#" ${pose}/full.gro ${pose}/correctBox.gro
sed -i ${pose}/solvated.gro -e "s/HOH/SOL/g"

echo "For ${pose} of ${PDBid}, the complex structure was solvated successfully."

done # poses

rm topol_amber.top

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, no pose could be solvated successfully. Deleting directory!"
cd ..
rm -r ${PDBid}
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
