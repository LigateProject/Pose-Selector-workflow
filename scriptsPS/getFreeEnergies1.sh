#!/bin/bash
set -e; shopt -s expand_aliases

module load gromacs/2023.2 gromacs=gmx

#PDBids=(6mj7 5yj8 5ovp 4mn3 5wxh 6j3p 4wn5 4zyf 6udu 6udt 5d3c 6dif 5dgw 6e9a)
PDBids=(5ovp 4mn3 5wxh 6j3p 4wn5 4zyf 6udu 6udt 5d3c 6dif 5dgw 6e9a)

PATH_TO_MDP=scriptsPS/mdp
solvationModel="pb"

for PDBid in ${PDBids[@]}; do
cd ${PDBid}

## update topology
ligandIdentifier=$(grep "Protein" topol_amber.top | tail -1 | cut -d " " -f 1)
### small organic molecule or any combinations or special cases
if grep -q "MOL " topol_amber.top
then
lineNumber=$(grep -n "MOL" topol_amber.top | tail -1 | cut -d ":" -f 1)
### peptide ligand
elif grep -q "${ligandIdentifier} " topol_amber.top
then
lineNumber=$(grep -n "${ligandIdentifier}" topol_amber.top | tail -1 | cut -d ":" -f 1)
fi
head -${lineNumber} topol_amber.top >> topol_rerun.top

for i in $(seq -w 01 08); do
cd run_${i}

# create .tpr file for complex only
## update .gro file
lineNumber=$(grep -n "SOL" *.gro | head -1 | cut -d ":" -f 1)
lineNumber=$((${lineNumber} - 1))
atomCount=$(head -2 simulation.gro | tail -1)
head -${lineNumber} simulation.gro >> rerunTPR.gro
tail -1 simulation.gro >> rerunTPR.gro
lineNumber=$((${lineNumber} - 2))
sed -i rerunTPR.gro -e "s/"${atomCount}"/"${lineNumber}"/g"

## create new index file
### system may be charged
gmx grompp -f ${PATH_TO_MDP}/dummy.mdp -c rerunTPR.gro -p ../topol_rerun.top -o selection.tpr -po selectionOut.mdp -maxwarn 1
### identify ligand correctly
#### small organic molecule or any combinations or special cases
val1=$(grep -n "; Compound" ../topol_rerun.top | cut -d ":" -f 1)
if grep -q "MOL " ../topol_rerun.top
then
val2=$(grep -n "MOL " ../topol_rerun.top | cut -d ":" -f 1)
#### peptide ligand
elif grep -q "${ligandIdentifier} " ../topol_rerun.top
then
val2=$(grep -n "${ligandIdentifier} " ../topol_rerun.top | cut -d ":" -f 1)
#### TODO add clauses for DNA and RNA ligands
fi
diff=$(( ${val2} - ${val1} ))

ligand='"ligand" molecule '${diff}''
ligand=\'${ligand}\'
protein='"protein" not molecule '${diff}''
protein=\'${protein}\'
complex='"complex" group System'
complex=\'${complex}\'

gmx_command="gmx select -s selection.tpr -on indexRerun.ndx -select ${ligand} ${protein} ${complex}"

echo "#!/bin/bash" >> createIndex.sh
echo "set -e" >> createIndex.sh
echo "" >> createIndex.sh
echo $gmx_command >> createIndex.sh
bash createIndex.sh
rm createIndex.sh

## use starting structure for .tpr instead
rm rerunTPR.gro
gmx trjconv -f simulation.xtc -s selection.tpr -n indexRerun.ndx -o rerunTPR.gro -dump 0 <<-eof
complex
eof
rm selection.tpr selectionOut.mdp

## grompp
gmx grompp -f ${PATH_TO_MDP}/simulationAnalysis.mdp -c rerunTPR.gro -p ../topol_rerun.top -n indexRerun.ndx -o rerun.tpr -po rerunOut.mdp -maxwarn 1

# prepare trajectory for binding free energy calculation
## select frames
firstFrame=100
lastFrame=100000
## convert frame numbers into simulated times
timeInterval=10.0 # in ps
startTime=$(echo ${firstFrame} ${timeInterval} | awk -v OFMT=%.5f '{print $1 * $2}')
endTime=$(echo ${lastFrame} ${timeInterval} | awk -v OFMT=%.5f '{print $1 * $2}')

## remove periodic boundary conditions from trajectory
gmx trjconv -f simulation.xtc -s rerun.tpr -n indexRerun.ndx -o ${solvationModel}_1.xtc -b ${startTime} -e ${endTime} -pbc nojump <<-eof
complex
complex
eof
gmx trjconv -f ${solvationModel}_1.xtc -s rerun.tpr -n indexRerun.ndx -o ${solvationModel}.xtc -pbc no -fit progressive <<-eof
complex
complex
eof
rm ${solvationModel}_1.xtc

cd ..
done
done
