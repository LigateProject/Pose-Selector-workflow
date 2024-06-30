#!/bin/bash
set -e; shopt -s expand_aliases

PDBids=(6mj7 5yj8 5ovp 4mn3 5wxh 6j3p 4wn5 4zyf 6udu 6udt 5d3c 6dif 5dgw 6e9a)

solvationModel="gb"
conversionFactor=4.184

for PDBid in ${PDBids[@]}; do
cd ${PDBid}
touch freeEnergySummary.txt

echo "Free energies as sum of gas and solvation entropy" >> freeEnergySummary.txt

for i in $(seq -w 01 08); do
cd run_${i}

# calculate binding free energy (more accurate)
freeEnergy=$(grep "TOTAL" ${solvationModel}_mmpbsa.dat | tail -1 | awk '{print $2;}')
freeEnergy=$(echo ${freeEnergy} ${conversionFactor} | awk -v OFMT=%.2f '{print $1 * $2}') # convert from kcal/mol to kJ/mol

standardDeviation=$(grep "TOTAL" ${solvationModel}_mmpbsa.dat | tail -1 | awk '{print $4;}')
standardDeviation=$(echo ${standardDeviation} ${conversionFactor} | awk -v OFMT=%.2f '{print $1 * $2}') # convert from kcal/mol to kJ/mol

echo "run_" ${i} ": " ${freeEnergy} " +- " ${standardDeviation} >> ../freeEnergySummary.txt

cd ..
done

echo "Free energies on the basis of the interaction entropy" >> freeEnergySummary.txt

for i in $(seq -w 01 08); do
cd run_${i}

# calculate binding free energy (interaction entropy)
freeEnergy=$(grep "G binding" ${solvationModel}_mmpbsa.dat | awk '{print $4;}')
freeEnergy=$(echo ${freeEnergy} ${conversionFactor} | awk -v OFMT=%.2f '{print $1 * $2}') # convert from kcal/mol to kJ/mol

standardDeviation=$(grep "G binding" ${solvationModel}_mmpbsa.dat | awk '{print $6;}')
standardDeviation=$(echo ${standardDeviation} ${conversionFactor} | awk -v OFMT=%.2f '{print $1 * $2}') # convert from kcal/mol to kJ/mol

echo "run_" ${i} ": " ${freeEnergy} " +- " ${standardDeviation} >> ../freeEnergySummary.txt

cd ..
done

for i in $(seq -w 01 08); do
cd run_${i}

## clean up
rm ${solvationModel}_filtered.xtc _GMXMMPBSA_* *.prmtop gmx_MMPBSA.log COMPACT_MMXSA_RESULTS.mmxsa ${solvationModel}_mmpbsa.csv
#rm indexRerun.ndx rerunTPR.gro rerunOut.mdp rerun.tpr  ${solvationModel}.* 
#rm indexRerun.ndx rerunTPR.gro rerunOut.mdp rerun.tpr  ${solvationModel}.* ${solvationModel}_filtered.xtc _GMXMMPBSA_* *.prmtop gmx_MMPBSA.log COMPACT_MMXSA_RESULTS.mmxsa ${solvationModel}_mmpbsa.csv

cd ..
done
done
