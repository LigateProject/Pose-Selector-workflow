#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

PATH_TO_MDP=scriptsPS/mdp
PATH_TO_SLURM=

PDBid=$(pwd | rev | cut -d "/" -f 1 | rev)

# the index file only needs to be created for every pose because the number of water molecules may differ
# TODO: adapt for loop such that ligand and ion identification are not repeated if it becomes a performance bottleneck.
for pose in pose_*; do
# create index file for output writting
gmx grompp -f ${PATH_TO_MDP}/dummy.mdp -c ${pose}/ions.gro -p topol_amber_${pose}.top -o ${pose}/selection.tpr -po ${pose}/selectionOut.mdp || true

## error catching
## selection.tpr must exist
if ! [ -f ${pose}/selection.tpr ]
then
echo "For ${pose} of ${PDBid}, the index file for output writing could not be generated (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

## identify ligand correctly
ligandIdentifier=$(grep "Protein" topol_amber_${pose}.top | tail -1 | cut -d " " -f 1)
### small organic molecule or any combinations or special cases
val1=$(grep -n "; Compound" topol_amber_${pose}.top | cut -d ":" -f 1)
if grep -q "MOL " topol_amber_${pose}.top
then
val2=$(grep -n "MOL " topol_amber_${pose}.top | cut -d ":" -f 1)
### peptide ligand
elif grep -q "${ligandIdentifier} " topol_amber_${pose}.top
then
val2=$(grep -n "${ligandIdentifier} " topol_amber_${pose}.top | cut -d ":" -f 1)
### TODO add clauses for DNA and RNA ligands
fi
diff=$(( ${val2} - ${val1} ))

## identify NA or CL ions that are part of the protein complex
stringProtein="group Protein or (group Ion and not (resname NA or resname CL))"
stringIons="(resname NA or resname CL)"
lineSolvent=$(grep -m1 -n "SOL" ${pose}/ions.gro | cut -d ":" -f 1)
### NA
linesNA=$(grep -n "NA      NA" ${pose}/ions.gro | cut -d ":" -f 1)
atomnrNA=$(grep "NA      NA" ${pose}/ions.gro | cut -c16-21)
for i in $(seq 1 $(echo $linesNA | wc -w)); do
if (( $(echo $linesNA | awk '{print $'${i}'}') < ${lineSolvent} ))
then
addition=$(echo ${atomnrNA} | awk '{print $'${i}'}')
stringProtein=${stringProtein}" or atomnr "${addition}
stringIons=${stringIons}" and not atomnr "${addition}
fi
done
### CL
linesCL=$(grep -n "CL      CL" ${pose}/ions.gro | cut -d ":" -f 1)
atomnrCL=$(grep "CL      CL" ${pose}/ions.gro | cut -c16-21)
for i in $(seq 1 $(echo $linesCL | wc -w)); do
if (( $(echo $linesCL | awk '{print $'${i}'}') < ${lineSolvent} ))
then
addition=$(echo ${atomnrCL} | awk '{print $'${i}'}')
stringProtein=${stringProtein}" or atomnr "${addition}
stringIons=${stringIons}" and not atomnr "${addition}
fi
done

## final formatting
if [[ ${stringIons: -1} != ")" ]]
then
stringIons="("${stringIons}")"
else
stringIons=${stringIons::-1}
stringIons=${stringIons:1}
fi
stringIons=" "${stringIons}
stringProtein="("${stringProtein}")"

## we have complex ions, but no other ions were added
noIons=0
if ! grep -q "NA      NA" ${pose}/ions.gro && ! grep -q "CL      CL" ${pose}/ions.gro
then
noIons=1
stringIons=""
stringProtein="(group Protein or group Ion)"
echo "For ${pose} of ${PDBid}, no ions were added previously."
elif grep -q "NA      NA" ${pose}/ions.gro && grep -q "CL      CL" ${pose}/ions.gro
then
if (( $(grep -n "NA      NA" ${pose}/ions.gro | tail -1 | cut -d ":" -f 1) < ${lineSolvent} )) &&  (( $(grep -n "CL      CL" ${pose}/ions.gro | tail -1 | cut -d ":" -f 1) < ${lineSolvent} ))
then
noIons=1
stringIons=""
stringProtein="(group Protein or group Ion)"
echo "For ${pose} of ${PDBid}, no ions were added previously."
fi
elif grep -q "NA      NA" ${pose}/ions.gro
then
if (( $(grep -n "NA      NA" ${pose}/ions.gro | tail -1 | cut -d ":" -f 1) < ${lineSolvent} ))
then
noIons=1
stringIons=""
stringProtein="(group Protein or group Ion)"
echo "For ${pose} of ${PDBid}, no ions were added previously."
fi
elif grep -q "CL      CL" ${pose}/ions.gro
then
if  (( $(grep -n "CL      CL" ${pose}/ions.gro | tail -1 | cut -d ":" -f 1) < ${lineSolvent} ))
then
noIons=1
stringIons=""
stringProtein="(group Protein or group Ion)"
echo "For ${pose} of ${PDBid}, no ions were added previously."
fi
fi

ligand='"ligand" molecule '${diff}''
ligand=\'${ligand}\'
protein='"protein" '${stringProtein}' and not molecule '${diff}'' # peptide ligands would be included in the group Protein
protein=\'${protein}\'
water='"water" resname SOL'
water=\'${water}\'
ions='"ions"'${stringIons}''
ions=\'${ions}\'
complex='"complex" molecule '${diff}' or '${stringProtein}''
complex=\'${complex}\'

# check if there are ions to look for
if ! [[ -z ${stringIons} ]]
then
stringIons=" or "${stringIons}
fi

waterAndIons='"water_and_ions" resname SOL'${stringIons}''
waterAndIons=\'${waterAndIons}\'

if (( ${noIons} == 1 ))
then
gmx_command="gmx select -s ${pose}/selection.tpr -on indexTraining_${pose}.ndx -select ${ligand} ${protein} ${water} ${complex} ${waterAndIons}"
else
gmx_command="gmx select -s ${pose}/selection.tpr -on indexTraining_${pose}.ndx -select ${ligand} ${protein} ${water} ${ions} ${complex} ${waterAndIons}"
fi

echo "#!/bin/bash" >> createIndex.sh
echo "set -e" >> createIndex.sh
echo "" >> createIndex.sh
echo $gmx_command >> createIndex.sh
bash createIndex.sh || true

## handle cases in which no ions are present (gmx select cannot match 'group "Ion"')
if ! [ -s indexTraining_${pose}.ndx ]
then
tailNumber=$(grep -n "gmx select," ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
lineNumber=$(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l)
tailNumber=$(echo ${lineNumber} ${tailNumber} | awk '{print $1-$2}')
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Cannot match 'group \"Ion\"'")
then
sed -i createIndex.sh -e "s/(group Protein or group Ion)/group Protein/g"
bash createIndex.sh || true
fi
fi

## error catching
## indexTraining.ndx must exist and must not be empty
if ! [ -s indexTraining_${pose}.ndx ]
then
echo "For ${pose} of ${PDBid}, the index file for output writing could not be generated (select error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

rm createIndex.sh ${pose}/selection.tpr ${pose}/selectionOut.mdp

# run grompp to get input .tpr file
gmx grompp -f ${PATH_TO_MDP}/simulation.mdp -c ${pose}/ions.gro -p topol_amber_${pose}.top -n indexTraining_${pose}.ndx -o ${pose}/simulation.tpr -po ${pose}/simulationOut.mdp || true
## error catching
## simulation.tpr must exist
if ! [ -f ${pose}/simulation.tpr ]
then
echo "For ${pose} of ${PDBid}, the TPR file could not be generated (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -r ${pose} topol_amber_${pose}.top
continue
fi

rm ${pose}/simulationOut.mdp

echo "For ${pose} of ${PDBid}, the TPR file has been generated successfully."

done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the TPR file could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
