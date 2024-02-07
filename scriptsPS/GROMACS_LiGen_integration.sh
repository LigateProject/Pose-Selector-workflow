#!/bin/bash
set -e

# required external software
module load python/3.11.1
module load forSTaGE/openbabel/3.1.1
module load forSTaGE/ambertools/21
module load forSTaGE/acpype/2022.7.21
module load gromacs/2023.2 gromacs=gmx

# STaGE itself
module load stage/1.0.0
export GMXLIB=$GMXDATA

INPUT_PATH=
PATH_TO_SCRIPTS=scriptsPS
FF=gaff2

PATH_TO_SLURM=

PDBid=$(pwd | rev | cut -d "/" -f 1 | rev)

echo "Generating topologies for ${PDBid} ..."

# pdb2gmx changes the atom numbering once again. Let's prepare for that.
if [ -s missingAtoms.txt ]
then
# convert atom numbers to unique descriptions
touch newMissingAtoms.txt
selection=""
numberOfFirstResidue=$(head -1 ${PDBid}.pdb | cut -c23-27)
searchTerms=$(tail -1 missingAtoms.txt)
for searchTerm in ${searchTerms}; do
findings=$(grep "${searchTerm}" ${PDBid}.pdb | cut -c$((12-${#searchTerm}))-12)
count=0
for finding in ${findings}; do
count=$((${count} + 1))
if [[ ${finding} == ${searchTerm} ]]
then
break
fi
done

chain=$(grep "${searchTerm}" ${PDBid}.pdb | head -${count} | tail -1 | cut -c22)
resnr=$(grep "${searchTerm}" ${PDBid}.pdb | head -${count} | tail -1 | cut -c23-27)
resnr=$(( ${resnr} - ${numberOfFirstResidue} + 1 ))
atomname=$(grep "${searchTerm}" ${PDBid}.pdb | head -${count} | tail -1 | cut -c14-17)
selection=$(echo ${selection}"(group Protein and chain "${chain}" and resnr "${resnr}" and name "${atomname}") or ")
done

selection=${selection::-4}

echo ${selection} >> newMissingAtoms.txt
mv newMissingAtoms.txt missingAtoms.txt
sed -i missingAtoms.txt -e "s/  )/)/g"
sed -i missingAtoms.txt -e "s/ )/)/g"
fi

file=${PDBid}.pdb

# ion residue naming might not meet pdb2gmx's expectations
# TODO: proper fix
sed -i ${file} -e "s/CA   CAS/CA    CA/g"

# check whether the PDB file contains ions or other substances that are not known to the force field used
unknownResidues=("MN" "NI" "CU" "CD" "FE" "CO" "SR")
unknownResidues2=("YCM")
unknownResidueFound=0

for unknownResidue in ${unknownResidues[@]}; do
if grep -q "${unknownResidue}    ${unknownResidue}" ${file} || grep -q "${unknownResidue}   ${unknownResidue}" ${file}
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, an unknown residue (${unknownResidue}) was found in the structure. Deleting folder as this structure should not be handled by the current workflow."
unknownResidueFound=1
break
fi
done

for unknownResidue in ${unknownResidues2[@]}; do
if grep -q "${unknownResidue}" ${file}
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, an unknown residue (${unknownResidue}) was found in the structure. Deleting folder as this structure should not be handled by the current workflow."
unknownResidueFound=1
break
fi
done

if (( ${unknownResidueFound} == 1 ))
then
exit 0
fi

# create protein topology
# AMBER99SB-ILDN
# TIP3P
## CONECT records between ions and amino acids are ignored!!!
# identify disulfide bonds
foundDisulfideBond=0
for i in $(seq 1 $(grep "SG" CONECT.txt | wc -l)); do
numberOfSGOccurences=$(grep "SG" CONECT.txt | head -${i} | tail -1 | awk -F"SG" '{print NF-1}')
positionOfFirstSG=$(grep "SG" CONECT.txt | head -${i} | tail -1 | awk -F"SG" '{print$1}')
positionOfFirstSG=${#positionOfFirstSG}
if (( ${numberOfSGOccurences} == 2 )) && (( ${positionOfFirstSG} < 10 ))
then
foundDisulfideBond=1
break
fi
done
if (( ${foundDisulfideBond} == 1 ))
then
gmx pdb2gmx -f $file -renum -ignh -chainsep ter -merge all <<-eof
6
1
eof

# catch long bonds
tailNumber=$(grep -n "Generating topologies for" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Long Bond")
then
echo "For ${PDBid}, the protein has long bonds. Deleting folder as this complex should not be handled by the current workflow."
cd ..
rm -r ${PDBid}
exit 0
fi

# merge all changes the residue numbering. Update missingResidues.txt and missingAtoms.txt accordingly.
# update residue numbering
firstResidues=$(grep "N3" topol.top | grep "N " | cut -c 18-25 | xargs)

## create lists of chains and the numbers of their initial residues
if $(grep -q "TER" ${PDBid}.pdb)
then
numberOfChains=$(grep -c "TER" ${PDBid}.pdb)
chainNames=""
for chain in $(seq 1 ${numberOfChains}); do
lineNumber=$(grep -n "TER" ${PDBid}.pdb | head -${chain} | tail -1 | cut -d ":" -f 1)
chainNames=$(echo ${chainNames}" "$(head -$((${lineNumber} - 1)) ${PDBid}.pdb | tail -1 | cut -c22))
done
else
numberOfChains=1
chainNames=$(grep -n "ATOM" ${PDBid}.pdb | head -1 | cut -d ":" -f 2 | cut -c22)
fi

## loop over files
filenames="missingResidues.txt missingAtoms.txt"
for filename in ${filenames}; do
selection=""
if [ -s ${filename} ]
then
missing=$(tail -1 ${filename})
if $(grep -q "TER" ${PDBid}.pdb)
then
chainName=""
residueNumber=""
## loop over file content
for word in ${missing}; do
## signal chain name
if [[ $(echo ${selection} | rev | cut -d " " -f 1 | rev) == "chain" ]]
then
chainName=${word}
elif [[ $(echo ${selection} | rev | cut -d " " -f 1 | rev) == "resnr" ]]
then
## formatting
if [[ ${filename} = "missingResidues.txt" ]]
then
word=$(echo ${word} | cut -d ")" -f 1)
fi
## calculate new residue number
residueNumber=$((${word} + $(echo ${firstResidues} | cut -d " " -f $(($(echo ${chainNames} | sed -e "s/ //g" | grep -b -o ${chainName} | cut -d ":" -f 1) + 1))) - 1))
## formatting
if [[ ${filename} = "missingResidues.txt" ]]
then
residueNumber=${residueNumber}")"
fi
word=${residueNumber}
## clear variables after calculating new residue number
chainName=""
residueNumber=""
fi
## append to file content
selection=${selection}" "${word}
done
else
selection=$(tail -1 ${filename})
fi
## remove useless chain identifiers from selection
for chainName in ${chainNames}; do
selection=$(echo ${selection} | sed -e "s/ and chain "${chainName}"//g")
done
## replace files
touch NEW_${filename}
echo ${selection} >> NEW_${filename}
mv NEW_${filename} ${filename}
fi
done

else
gmx pdb2gmx -f $file -renum -ignh -chainsep ter <<-eof
6
1
eof
fi

# catch long bonds
tailNumber=$(grep -n "Generating topologies for" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Long Bond")
then
echo "For ${PDBid}, the protein has long bonds. Deleting folder as this complex should not be handled by the current workflow."
cd ..
rm -r ${PDBid}
exit 0
fi

# generate topology and .gro files for ligand
## only the cleaned .mol2 file has the correct atom naming
cp ${INPUT_PATH}/${PDBid}/${PDBid}_ligand_cleaned.mol2 ligand.mol2
obabel "ligand.mol2" -O "${PDBid}.mol2"

## correct residue names in cleaned .mol2 file
residueNames1=$(awk '/@<TRIPOS>ATOM/{flag=1; next}/@<TRIPOS>/{flag=0} flag' ${PDBid}.mol2 | awk '{print $8}')
residueNames1=(${residueNames1})
for residueName1 in ${residueNames1[@]}; do
if ! $(echo "${residueName1}" | grep -q "\*")
then
if (( ${#residueName1} > 3 ))
then
residueName1length=${#residueName1}
spacesToAdd=$(( ${residueName1length} - 3 ))
lineOfSpaces=$(printf ' %.0s' $(seq 1 $spacesToAdd))
residueName1New="$(echo ${residueName1} | cut -c1-3)${lineOfSpaces}"
sed -i ${PDBid}.mol2 -e "s/${residueName1}/${residueName1New}/g"
fi
fi
done

## extract poses
numberOfPoses=$(grep -c "<TRIPOS>MOLECULE" ${INPUT_PATH}/${PDBid}/ligen_poses.mol2)
numberOfPoses=$(( ${numberOfPoses} - 1 )) # counting from 0

for i in $(seq 0 ${numberOfPoses}); do
mkdir pose_${i}
python3 ${PATH_TO_SCRIPTS}/extractPose.py <<-eof
${INPUT_PATH}/${PDBid}/ligen_poses.mol2
${i}
pose_${i}/ligand.mol2
eof
obabel "pose_${i}/ligand.mol2" -O "pose_${i}/${PDBid}.mol2"
rm pose_${i}/ligand.mol2
done

## correct atom and residue naming for poses
startLineRef=$(grep -n "@<TRIPOS>ATOM" ${PDBid}.mol2 | cut -d ":" -f 1)
endLineRef=$(grep -Rn "@<TRIPOS>" ${PDBid}.mol2 | awk 'c&&!--c;/>ATOM/{c=1}' | cut -d ":" -f 1)

for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

touch pose_${i}/correctPose.mol2
startLinePose=$(grep -n "@<TRIPOS>ATOM" pose_${i}/${PDBid}.mol2 | cut -d ":" -f 1)
endLinePose=$(grep -Rn "@<TRIPOS>" pose_${i}/${PDBid}.mol2 | awk 'c&&!--c;/>ATOM/{c=1}' | cut -d ":" -f 1)
head -${startLinePose} pose_${i}/${PDBid}.mol2 >> pose_${i}/correctPose.mol2
for j in $(seq 1 $(( ${endLineRef} - ${startLineRef} - 1 ))); do
fileContentRef=$(head -$(( ${startLineRef} + ${j} )) ${PDBid}.mol2 | tail -1)
fileContentPose=$(head -$(( ${startLinePose} + ${j} )) pose_${i}/${PDBid}.mol2 | tail -1)
residueNameRef=$(echo "${fileContentRef}" | cut -c59-63 | sed -e "s/\*/\\\*/g")
residueNamePose=$(echo "${fileContentPose}" | cut -c59-63 | sed -e "s/\*/\\\*/g")
atomNameRef=$(echo "${fileContentRef}" | cut -c9-11)
atomNamePose=$(echo "${fileContentPose}" | cut -c9-11)
echo "${fileContentPose}" | sed -e "s/${residueNamePose}/${residueNameRef}/g" | sed -e "s/${atomNamePose}/${atomNameRef}/g" >> pose_${i}/correctPose.mol2
done
tail +${endLinePose} pose_${i}/${PDBid}.mol2 >> pose_${i}/correctPose.mol2
mv pose_${i}/correctPose.mol2 pose_${i}/${PDBid}.mol2
done
echo "For ${PDBid}, the atom and residue naming of the individual poses has been corrected."

# check ligand
## initialise
count=0
countPeptide=0
countDNA=0
countRNA=0
residueCount=0
previousResidueName=""
isPeptideLigand=0
## determine number of atoms
patternNumber=$(grep "TRIPOS" ${PDBid}.mol2 | grep -n "@<TRIPOS>ATOM" | cut -d ":" -f 1)
patternNumber=$(( ${patternNumber} + 1 ))
pattern=$(grep -m${patternNumber} "TRIPOS" ${PDBid}.mol2 | tail -1)
numberOfAtoms=$(awk '/@<TRIPOS>ATOM$/ {b=NR; next} /'${pattern}'$/ {print NR-b-1; exit}' ${PDBid}.mol2)
## determine range for loop over atoms
lineNumber=$(grep -n "@<TRIPOS>ATOM" ${PDBid}.mol2 | cut -d ":" -f 1)
lineNumber2=$(( ${lineNumber} + ${numberOfAtoms} ))
## loop over atoms
for i in $(seq -w $(( ${lineNumber} + 1 )) ${lineNumber2}); do
content=$(head -n ${i} ${PDBid}.mol2 | tail -n 1 | cut -c59-63)
# remove spaces
if ! $(echo "${content}" | grep -q "\*\*\*\*")
then
content=$(echo ${content})
fi
### check number of residues in ligand
### blind spot: a chain of the same type of residue is considered one residue
### however, pdb2gmx has the same blind spot and will refuse to parameterise a chain containing only one type of residue
if (( ${residueCount} < 2 ))
then
### in the dataset, all hydrogen atoms have a residue name containing ****
### therefore, residue names containing **** must not increase the residue count
if [[ ${content} != ${previousResidueName} ]]  && ! $(echo "${content}" | grep -q "\*\*\*\*")
then
residueCount=$(( ${residueCount} + 1 ))
fi
previousResidueName=${content}
fi
### check whether all residues are standard amino acids or nucleobases
### search for the residue name with a space at the end to avoid issues like finding residue "NVA" due to the existence of "NVAL"
### exclude residue names containing numbers
contentCheck=$(echo ${content} | tr -d -c 0-9)
if grep -q "${content} " ${GMXDATA}/amber99sb-ildn.ff/aminoacids.rtp && (( ${#contentCheck} == 0 ))
then
count=$(( ${count} + 1 ))
countPeptide=$(( ${countPeptide} + 1 ))
### ugly fix: some small organic compound is abbreviated the same way as an DNA nucleotide
### TODO: find a proper fix
elif grep -q "${content} " ${GMXDATA}/amber99sb-ildn.ff/dna.rtp && [[ ${content} != "DAN" ]] && (( ${#contentCheck} == 0 ))
then
count=$(( ${count} + 1 ))
countDNA=$(( ${countDNA} + 1 ))
### ugly fix: some small organic compound is abbreviated the same way as an RNA nucleotide
### TODO: find a proper fix
elif grep -q "${content} " ${GMXDATA}/amber99sb-ildn.ff/rna.rtp && [[ ${content} != "RUN" ]] && (( ${#contentCheck} == 0 ))
then
count=$(( ${count} + 1 ))
countRNA=$(( ${countRNA} + 1 ))
### in the PS dataset, all hydrogen atoms have **** as residue name => do not consider when identifying molecule type
elif $(echo "${content}" | grep -q "\*\*\*\*")
then
#### if the heavy atoms all belong to one kind of molecule, the group of hydrogens at the end of the file will not change the classification
if (( ${count} == ${countPeptide} ))
then
countPeptide=$(( ${countPeptide} + 1 ))
elif (( ${count} == ${countDNA} ))
then
countDNA=$(( ${countDNA} + 1 ))
elif (( ${count} == ${countRNA} ))
then
countRNA=$(( ${countRNA} + 1 ))
fi
#### we later on check that count equals the number of atoms in the system => increase it for hydrogens, too
count=$(( ${count} + 1 ))
else
break
fi
#### fall back to GAFF if we have a Frankenstein ligand combining amino acids, DNA and/or RNA
if (( ${count} > ${countPeptide} && ${countPeptide} > 0 )) || (( ${count} > ${countDNA} && ${countDNA} > 0 )) ||  (( ${count} > ${countRNA} && ${countRNA} > 0 ))
then
echo "${PDBid}: Mixture of peptide, DNA and/or RNA. Falling back to GAFF2 because individual specialised force fields cannot be combined easily."
break
fi
done
## check additional limitations
### the AMBER force field shipped with GROMACS has no entries for zwitterionic isolated/single amino acids => fall back to GAFF2
if (( ${residueCount} < 2 && ${countPeptide} > 0 ))
then
countPeptide=$(( ${countPeptide} - 1 ))
echo "${PDBid}: Only one residue. Falling back to GAFF2 because the protein force field does not contain parameters for zwitterionic amino acids."
fi

## rename fluorine atoms if unconventional naming is used, e.g. FAA, FAB should be changed to F01, F02
for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

if grep -qE -e '^.{8}F' pose_${i}/${PDBid}.mol2
then
countFluorine=1
while IFS= read -r line; do
    if [[ ${line:8:5} =~ F.... ]]; then
        replacement="F$(printf "%-4d" $countFluorine)"
        line="${line:0:8}$replacement${line:13}"
        countFluorine=$((countFluorine+1))
    fi
    echo "$line"
done < pose_${i}/${PDBid}.mol2 > pose_${i}/tmp_mol2
cat pose_${i}/tmp_mol2 > pose_${i}/${PDBid}.mol2
rm pose_${i}/tmp_mol2
fi
done

# treatment of peptide, DNA or RNA ligand
if (( ${countPeptide} == ${numberOfAtoms} ))
then
isPeptideLigand=1
echo "${PDBid}: PEPTIDE, DNA or RNA as LIGAND!"
for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

obabel "pose_${i}/${PDBid}.mol2" -O "pose_${i}/${PDBid}_ligand.pdb"

# sometimes, OpenBabel uses atom names that are incompatible with GROMACS
sed -i pose_${i}/${PDBid}_ligand.pdb -e "s/O1 /OXT/g"
sed -i pose_${i}/${PDBid}_ligand.pdb -e "s/O30/OXT/g"
sed -i pose_${i}/${PDBid}_ligand.pdb -e "s/O48/OXT/g"
sed -i pose_${i}/${PDBid}_ligand.pdb -e "s/O65/OXT/g"
sed -i pose_${i}/${PDBid}_ligand.pdb -e "s/OE /OXT/g"

# not happy about using -ignh
set +e
gmx pdb2gmx -f pose_${i}/${PDBid}_ligand.pdb -renum -ignh -p pose_${i}/ligand.top -i pose_${i}/posre_Ligand.itp -o pose_${i}/ligand.gro <<-eof
6
1
eof
set -e

# catch potential pdb2gmx errors
if ! [ -s pose_${i}/ligand.top ]
then
echo "For pose_${i} of ${PDBid}, pdb2gmx could not generate a topology for the peptide, DNA or RNA ligand. Deleting folder as this complex should not be handled by the current workflow."
rm -r pose_${i}
continue
fi

# catch long bonds
tailNumber=$(grep -n "pdb2gmx" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Long Bond")
then
echo "For pose_${i} of ${PDBid}, the peptide, DNA or RNA ligand has long bonds. Deleting folder as this complex should not be handled by the current workflow."
rm -r pose_${i}
continue
fi
done

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the ligand topology could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
exit 0
fi

# print .gro file for protein-ligand complex
for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

python3 $PATH_TO_SCRIPTS/printComplexGroFile.py <<-eof
conf.gro
pose_${i}/ligand.gro
pose_${i}/full.gro
eof

gmx editconf -f pose_${i}/full.gro -o pose_${i}/full.pdb

# add TER signal at the end of all protein chains
OC2=$(grep -c "OC2" pose_${i}/full.pdb)

for h in $(seq 1 ${OC2}); do
lineNumber=$(grep -n "OC2" pose_${i}/full.pdb | head -${h} | tail -1 | cut -d ":" -f 1)
lineNumber2=$(wc -l pose_${i}/full.pdb | cut -d " " -f 1)
head -${lineNumber} pose_${i}/full.pdb >> pose_${i}/full_final.pdb
echo "TER" >> pose_${i}/full_final.pdb
tail -$(( ${lineNumber2} - ${lineNumber} )) pose_${i}/full.pdb >> pose_${i}/full_final.pdb
mv pose_${i}/full_final.pdb pose_${i}/full.pdb
done

# add TER after small organic ligand
## peptide ligand already has a TER statement
if (( ${isPeptideLigand} == 0 ))
then
if $(grep -q "HOH" pose_${i}/full.pdb)
then
lineNumber=$(grep -n "HOH" pose_${i}/full.pdb | head -1 | cut -d ":" -f 1)
lineNumber=$((${lineNumber} - 1))
if ! $(head -${lineNumber} pose_${i}/full.pdb | tail -1 | grep -q "TER")
then
lineNumber2=$(wc -l pose_${i}/full.pdb | cut -d " " -f 1)
head -${lineNumber} pose_${i}/full.pdb >> pose_${i}/full_final.pdb
echo "TER" >> pose_${i}/full_final.pdb
tail -$(( ${lineNumber2} - ${lineNumber} )) pose_${i}/full.pdb >> pose_${i}/full_final.pdb
mv pose_${i}/full_final.pdb pose_${i}/full.pdb
fi
## We may already have a TER statement because it's the end of the PDB file.
## However, duplicated TER statements are harmless in this conditional branch.
else
echo "TER" >> pose_${i}/full.pdb
fi
fi

# clean up
rm pose_${i}/${PDBid}_ligand.pdb pose_${i}/${PDBid}.mol2 pose_${i}/ligand.gro pose_${i}/posre*itp pose_${i}/full.gro

# pdb2gmx has trouble identifying ions if their residue names start with numbers
# only observed for one structure so far
# TODO: proper fix if observed more frequently
sed -i pose_${i}/full.pdb -e "s/10ZN/  ZN/g"
sed -i pose_${i}/full.pdb -e "s/11ZN/  ZN/g"
sed -i pose_${i}/full.pdb -e "s/12ZN/  ZN/g"
sed -i pose_${i}/full.pdb -e "s/0ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/1ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/2ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/3ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/4ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/5ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/6ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/7ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/8ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/9ZN/ ZN/g"
sed -i pose_${i}/full.pdb -e "s/8CA/ CA/g"
sed -i pose_${i}/full.pdb -e "s/4NA/ NA/g"
sed -i pose_${i}/full.pdb -e "s/6NA/ NA/g"
sed -i pose_${i}/full.pdb -e "s/3K / K /g"
sed -i pose_${i}/full.pdb -e "s/4K / K /g"
sed -i pose_${i}/full.pdb -e "s/2MG/ MG/g"
sed -i pose_${i}/full.pdb -e "s/3MG/ MG/g"

# final topology
## CONECT records between ions and amino acids are ignored!!!
## As we do not let pdb2gmx regenerate H, the SS-bond pattern should be preserved.
if (( ${foundDisulfideBond} == 1 ))
then
### if the ligand is a peptide, we do not want it to be merged with the rest of the protein
### it is safe to assume that the ligand is always the last protein chain
### merge interactively and answer yes until you are asked about the ligand chain; then say no
if (( ${isPeptideLigand} == 1 ))
then
numberOfChainsToBeMerged=$(grep -c "TER" pose_${i}/full.pdb)
numberOfChainsToBeMerged=$(( ${numberOfChainsToBeMerged} - 3 ))
touch createTopology.sh
echo "gmx pdb2gmx -f pose_${i}/full.pdb -p pose_${i}/topol_amber.top -i pose_${i}/posre.itp -o pose_${i}/full.gro -chainsep ter -merge interactive <<-eof" >> createTopology.sh
echo "6" >> createTopology.sh
echo "1" >> createTopology.sh
for c in $(seq 1 ${numberOfChainsToBeMerged}); do
echo "y" >> createTopology.sh
done
echo "n" >> createTopology.sh
echo "eof" >> createTopology.sh
### For reasons I don't understand, pdb2gmx can end up using the wrong type of histidine (e.g. HISD vs. HISE)
### even if hydrogen atoms are present in the PDB file.
### For now, we exclude poses for which pdb2gmx wants to use a different type of HIS.
### TODO: read out first HIS pattern and enforce it in subsequent poses or code a GROMACS fix (preferred)
set +e
bash createTopology.sh
set -e
rm createTopology.sh
### standard treatment of proteins with disulfide bonds and non-peptide ligands
else
set +e
gmx pdb2gmx -f pose_${i}/full.pdb -p pose_${i}/topol_amber.top -i pose_${i}/posre.itp -o pose_${i}/full.gro -chainsep ter -merge all <<-eof
6
1
eof
set -e
fi
else
set +e
gmx pdb2gmx -f pose_${i}/full.pdb -p pose_${i}/topol_amber.top -i pose_${i}/posre.itp -o pose_${i}/full.gro -chainsep ter <<-eof
6
1
eof
set -e
fi

# The topology for the pose is only written if pdb2gmx succeeds.
if ! [ -s pose_${i}/topol_amber.top ]
then
echo "For pose_${i} of ${PDBid}, pdb2gmx failed to generate a topology. Deleting folder as this complex should not be handled by the current workflow."
rm -r pose_${i}
continue
fi

existingPose="0"

for h in $(seq 0 ${numberOfPoses}); do
if [ -d pose_${h} ]
then
existingPose=${h}
break
fi
done

# this assumes that pdb2gmx uses the same protonation pattern for all poses (dangerous with -ignh)
if (( ${i} > ${existingPose} ))
then
# check if topology differs from top level topology
# TODO: go with individual topologies for each pose if we have a peptide ligand (only if the test below fails very often)
# TODO: can we make this check more rigorous?
## compare
difference=$(diff pose_${existingPose}/ligand.top pose_${i}/ligand.top | wc -l)
if (( ${difference} > 20 ))
then
echo "For pose_${i} of ${PDBid}, the topology created by pdb2gmx differs from the topology of pose ${existingPose}. Deleting folder as this complex should not be handled by the current workflow."
rm -r pose_${i}
continue
fi
rm pose_${i}/*.top
rm pose_${i}/*.itp
fi

# final clean up
rm pose_${i}/full.pdb
done

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the ligand topology could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
exit 0
fi

rm *.mol2 conf.gro topol* posre*itp
mv pose_${existingPose}/*.top .
mv pose_${existingPose}/*.itp .
sed -i *.top -e "s=pose_${existingPose}/==g"
sed -i *.itp -e "s=pose_${existingPose}/==g"
rm ligand.top

# treatment of general ligand
else
echo "${PDBid}: Assuming a SMALL ORGANIC LIGAND!"
stage=$(which stage.py)
python3 $stage -i "${PDBid}.mol2" -o "${PDBid}_stage" --forcefields $FF

# error catching
## if STaGE fails, try again with different total charge
lineNumberForCheckingStage=$(grep -n "Assuming a SMALL ORGANIC LIGAND!" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
lineNumberForCheckingStage=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${lineNumberForCheckingStage} ))
if $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Error running generator for gaff2") || $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Error removing OPLS directory")
then
rm -r *${PDBid}_stage*
python3 $stage -i "${PDBid}.mol2" -o "${PDBid}_stage" --forcefields $FF --alternativeChargeRounding
fi
# if it still fails, give up
lineNumberForCheckingStage=$(grep -n "Assuming a SMALL ORGANIC LIGAND!" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
lineNumberForCheckingStage=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${lineNumberForCheckingStage} ))
if (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error running generator for gaff2") == 2 )) || (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error removing OPLS directory") == 2 ))
then
echo "For ${PDBid}, STaGE failed to generate a ligand topology. Deleting folder as this ligand should not be handled by the current workflow."
cd ..
rm -r ${PDBid}
exit 0
elif (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error running generator for gaff2") == 1 )) || (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error removing OPLS directory") == 1 ))
then
echo "For ${PDBid}, STaGE worked with the alternative total charge."
else
echo "For ${PDBid}, STaGE successfully generated the ligand topology with the total charge from the input MOL2 file."
fi

# sort files
sed -i ${PDBid}_stage_${FF}/${PDBid}_stage.itp -e "s/${PDBid}_stage/ligand/g"
mv ${PDBid}_stage.gro ligand.gro
mv posre_${PDBid}_stage.itp posre_Ligand.itp
mv ${PDBid}_stage_${FF}/${PDBid}_stage.itp ligand_premature.itp

python3 $PATH_TO_SCRIPTS/topologySplitter.py <<-eof
ligand_premature.itp
eof

# write topology summaries to be used for simulations
python3 $PATH_TO_SCRIPTS/writeTopologySummarySingleLigand.py

# clean-up
rm -r ${PDBid}_stage_${FF} 
rm *.mol2 ligand_premature.itp topol.top ligand.gro
rm topol_ligandInWater.top

# create .gro files for poses and clean up
for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

obabel "pose_${i}/${PDBid}.mol2" -O "pose_${i}/ligand.gro"
rm pose_${i}/*.mol2
# print .gro file for protein-ligand complex
python3 $PATH_TO_SCRIPTS/printComplexGroFile.py <<-eof
conf.gro
pose_${i}/ligand.gro
pose_${i}/full.gro
eof

# clean-up
rm pose_${i}/ligand.gro
done

# clean-up
rm conf.gro

fi

# check whether residues very close to the ligand had to be modelled
if [ -s missingResidues.txt ]
then
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
## create .tpr file to produce index file (gmx select requires a .tpr file containing velocities)
gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c pose_${i}/full.gro -r pose_${i}/full.gro -p topol_amber.top -o pose_${i}/check.tpr -po pose_${i}/checkout.mdp -maxwarn 2 || true

### catch potential grompp segfault
if $(tail -1 ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "core dumped")
then
echo "For pose_${i} of ${PDBid}, executing gmx grompp resulted in a segmentation fault/core dump. Deleting folder as this complex should not be handled by the current workflow."
rm -r pose_${i}/
continue
fi

### catch other grompp errors
if [ ! -s pose_${i}/check.tpr ]
then
echo "For pose_${i} of ${PDBid}, the TPR file could not be generated (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -r pose_${i}/
continue
fi

fi
done

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the TPR file could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
exit 0
fi

## produce index file (gmx select errors out on composed command => print it to different bash script first)
ligandIdentifier=$(grep "Protein" topol_amber.top | tail -1 | cut -d " " -f 1)
val1=$(grep -n "; Compound" topol_amber.top | cut -d ":" -f 1)
if grep -q "MOL " topol_amber.top
then
val2=$(grep -n "MOL " topol_amber.top | cut -d ":" -f 1)
elif grep -q "${ligandIdentifier} " topol_amber.top
then
val2=$(grep -n "${ligandIdentifier} " topol_amber.top | cut -d ":" -f 1)
fi
diff=$(( ${val2} - ${val1} ))
ligand='"ligand" molecule '${diff}''
ligand=\'${ligand}\'

neighbours=$(tail -1 missingResidues.txt)
neighbours='"modelledResidues" '${neighbours}'' # promod3 doesn't model H
# keyword chain does not work correctly
if grep -q "Protein_chain" topol_amber.top
then
numberOfChains=$(grep -c "Protein_chain" topol_amber.top)
numberOfChains=$(( ${numberOfChains}/2 ))
for c in $(seq 1 ${numberOfChains}); do
neighbours=$(echo ${neighbours} | sed -e "s/chain $(grep "Protein_chain" topol_amber.top | head -${c} | tail -1 | cut -d "r" -f 2 | cut -d "_" -f 3 | cut -c1)/molecule ${c}/g")
done
else
neighbours=$(echo ${neighbours} | sed -e "s/chain . and //g")
fi
neighbours=\'${neighbours}\'

existingPose="0"

for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i}/check.tpr ]
then
existingPose=${i}
break
fi
done

gmx_command="gmx select -s pose_${existingPose}/check.tpr -on check.ndx -select ${ligand} ${neighbours}"

echo "#!/bin/bash" >> createIndex.sh
echo "set -e" >> createIndex.sh
echo "" >> createIndex.sh
echo $gmx_command >> createIndex.sh
bash createIndex.sh
rm createIndex.sh

## calculate minimum distance between ligand and modelled residues
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
gmx mindist -f pose_${i}/full.gro -s pose_${i}/check.tpr -n check.ndx -od pose_${i}/check.xvg <<-eof
ligand
modelledResidues
eof
distance=$(tail -1 pose_${i}/check.xvg | rev | cut -d " " -f 1 | rev) # distance in nm
distance=$(echo ${distance} 10 | awk '{print $1 * $2}') # distance in A
echo "For pose_${i} of ${PDBid}, the minimum distance between the ligand and modelled residues (in A) amounts to: ${distance}"
if (( $(echo "${distance} < 3" | bc -l) ))
then
echo "For pose_${i} of ${PDBid}, modelled residues are very close to the bound ligand (less than 3 A). Deleting folder as this structure should not be handled by the current workflow!"
rm -r pose_${i}
continue
fi

## clean up
rm pose_${i}/checkout.mdp pose_${i}/check.tpr pose_${i}/check.xvg
fi
done
rm check.ndx

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, modelled residues are very close to the bound ligand (less than 3 A) for all poses. Deleting directory!"
cd ..
rm -r ${PDBid}
exit 0
fi

else
echo "For ${PDBid}, no residues needed to be modelled."
fi


# check whether residues very close to the ligand had to be modelled
if [ -s missingAtoms.txt ]
then
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
## create .tpr file to produce index file (gmx select requires a .tpr file containing velocities)
gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c pose_${i}/full.gro -r pose_${i}/full.gro -p topol_amber.top -o pose_${i}/check.tpr -po pose_${i}/checkout.mdp -maxwarn 2 || true

### catch potential grompp segfault
if $(tail -1 ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "core dumped")
then
echo "For pose_${i} of ${PDBid}, executing gmx grompp resulted in a segmentation fault/core dump. Deleting folder as this complex should not be handled by the current workflow."
rm -r pose_${i}/
continue
fi

### catch other grompp errors
if [ ! -s pose_${i}/check.tpr ]
then
echo "For pose_${i} of ${PDBid}, the TPR file could not be generated (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -r pose_${i}/
continue
fi

fi
done

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, the TPR file could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -r ${PDBid}
exit 0
fi

## produce index file (gmx select errors out on composed command => print it to different bash script first)
ligandIdentifier=$(grep "Protein" topol_amber.top | tail -1 | cut -d " " -f 1)
val1=$(grep -n "; Compound" topol_amber.top | cut -d ":" -f 1)
if grep -q "MOL " topol_amber.top
then
val2=$(grep -n "MOL " topol_amber.top | cut -d ":" -f 1)
elif grep -q "${ligandIdentifier} " topol_amber.top
then
val2=$(grep -n "${ligandIdentifier} " topol_amber.top | cut -d ":" -f 1)
fi
diff=$(( ${val2} - ${val1} ))
ligand='"ligand" molecule '${diff}''
ligand=\'${ligand}\'

neighbours=$(tail -1 missingAtoms.txt)
neighbours='"modelledAtoms" '${neighbours}'' # promod3 doesn't model H
# keyword chain does not work correctly
if grep -q "Protein_chain" topol_amber.top
then
numberOfChains=$(grep -c "Protein_chain" topol_amber.top)
numberOfChains=$(( ${numberOfChains}/2 ))
for c in $(seq 1 ${numberOfChains}); do
neighbours=$(echo ${neighbours} | sed -e "s/chain $(grep "Protein_chain" topol_amber.top | head -${c} | tail -1 | cut -d "P" -f 2 | cut -d "_" -f 3 | cut -c1)/molecule ${c}/g")
done
else
neighbours=$(echo ${neighbours} | sed -e "s/chain . and //g")
fi
neighbours=\'${neighbours}\'

existingPose="0"

for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i}/check.tpr ]
then
existingPose=${i}
break
fi
done

gmx_command="gmx select -s pose_${existingPose}/check.tpr -on check.ndx -select ${ligand} ${neighbours}"

echo "#!/bin/bash" >> createIndex.sh
echo "set -e" >> createIndex.sh
echo "" >> createIndex.sh
echo $gmx_command >> createIndex.sh
bash createIndex.sh
rm createIndex.sh

## calculate minimum distance between ligand and modelled residues
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
gmx mindist -f pose_${i}/full.gro -s pose_${i}/check.tpr -n check.ndx -od pose_${i}/check.xvg <<-eof
ligand
modelledAtoms
eof
distance=$(tail -1 pose_${i}/check.xvg | rev | cut -d " " -f 1 | rev) # distance in nm
distance=$(echo ${distance} 10 | awk '{print $1 * $2}') # distance in A
echo "For pose_${i} of ${PDBid}, the minimum distance between the ligand and modelled atoms (in A) amounts to: ${distance}"
if (( $(echo "${distance} < 3" | bc -l) ))
then
echo "For pose_${i} of ${PDBid}, modelled atoms are very close to the bound ligand (less than 3 A). Deleting folder as this protein should not be handled by the current workflow!"
rm -r pose_${i}
continue
fi

## clean up
rm pose_${i}/checkout.mdp pose_${i}/check.tpr pose_${i}/check.xvg
fi
done
rm check.ndx

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${PDBid}, modelled atoms are very close to the bound ligand (less than 3 A) for all poses. Deleting directory!"
cd ..
rm -r ${PDBid}
exit 0
fi

else
echo "For ${PDBid}, no atoms in residues needed to be modelled."
fi


# final clean up
rm ${PDBid}.pdb missingResidues.txt missingAtoms.txt

# check if CONECT records have been ignored (metals are always ignored; only check disulfide bonds)
if [ -s CONECT.txt ]
then

## identify bonds given by CONECT statements
numberOfBonds=$(grep "bond" CONECT.txt | wc -l)
numberOfDisulfideBonds=0
numberOfDisulfideBondsTopology=0
for b in $(seq 1 ${numberOfBonds}); do
bond=$(grep -m${b} "bond" CONECT.txt | tail -1)
atom1=$(echo ${bond} | cut -d " " -f 2)
atom2=$(echo ${bond} | cut -d " " -f 3)

## identify type of atom
### we always have to take the atom type that is listed directly after the first occurrence of the respective atom number in CONECT.txt
element=$(( $(grep -m1 "${atom1}" CONECT.txt | wc -w) - $(grep -m1 "${atom1}" CONECT.txt | awk -F "${atom1}" '{print $2}' | wc -w) + 1 ))
type1=$(grep -m1 "${atom1}" CONECT.txt  | awk '{print $'${element}'}')
element=$(( $(grep -m1 "${atom2}" CONECT.txt | wc -w) - $(grep -m1 "${atom2}" CONECT.txt | awk -F "${atom2}" '{print $2}' | wc -w) + 1 ))
type2=$(grep -m1 "${atom2}" CONECT.txt  | awk '{print $'${element}'}')

## count number of disulfide bonds
if [[ ${type1} == "SG" && ${type2} == "SG" ]]
then
numberOfDisulfideBonds=$((${numberOfDisulfideBonds} + 1))
fi
done # bonds

## identify file containing protein topology
if $(ls | grep -q "_Protein")
then
files=$(ls topol_*.itp)
else
files="topol_amber.top"
fi

for file in ${files}; do
## identify all SG in CYS
allSG=""
for s in $(seq 1 $(grep "CYS     SG" ${file} | wc -l)); do
allSG=${allSG}" "$(grep "CYS     SG" ${file} | head -${s} | tail -1 | cut -c1-6 | tr -c -d 0-9)
done

## identify lines where bonds are listed
lineNumber1=$(grep -n "\[ bonds \]" ${file} | cut -d ":" -f 1)
patternNumber=$(grep "\[" ${file} | grep -n "\[ bonds \]" | cut -d ":" -f 1)
patternNumber=$(( ${patternNumber} + 1 ))
pattern=$(grep -m${patternNumber} "\[" ${file} | tail -1 | sed -e "s/\[//g" | sed -e "s/\]//g")
lineNumber2=$(grep -n "\[${pattern}\]" ${file} | cut -d ":" -f 1)

## extract bonds from topology
bondsTopology=$(head -$((${lineNumber2} - 2)) ${file} | tail -$((${lineNumber2} - ${lineNumber1} - 3)))

## search for disulfide bonds
for s in ${allSG}; do
for d in ${allSG}; do
### we have to exclude bonds that have a different first atom whose number just happens to end on the same digits as the atom number of a SG in a disulfide bond
### SG can never be the first atom in a topology such that we can add a leading space to fulfill the above requirement
### similarly, we have to exclude bonds that have a different second atom whose number just happens to begin with the same digits as the atom number of a SG in a disulfide bond
### the bond type is always indicated such that we can add a space after the second atom
if $(echo ${bondsTopology} | grep -q " ${s} ${d} ") && (( ${s} != ${d} ))
then
numberOfDisulfideBondsTopology=$((${numberOfDisulfideBondsTopology} + 1))
fi
done
done

done # files

## PDB CONECT statements list all disulfide bonds twice
numberOfDisulfideBondsTopology=$((${numberOfDisulfideBondsTopology} + ${numberOfDisulfideBondsTopology}))

if (( ${numberOfDisulfideBonds} != ${numberOfDisulfideBondsTopology} ))
then
echo "For ${PDBid}, the number of disulfide bonds in the PDB CONECT statements and in the GROMACS topology is not identical. Check topologies manually before continuing!"
touch incorrect_number_of_disulfide_bond
fi

fi

rm CONECT.txt residuesNotRemodelled.txt

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload forSTaGE/openbabel/3.1.1
module unload forSTaGE/ambertools/21
module unload forSTaGE/acpype/2022.7.21
module unload gromacs/2023.2

export GMXLIB=""
