#!/bin/bash
set -e

# make sure we start from a clean environment (mainly needed if the previous complex was deleted)
modules=(biopython/1.81 boost_own/1.82.0 forSTaGE/acpype/2022.7.21 forSTaGE/ambertools/21 forSTaGE/openbabel/3.1.1 gromacs/2023.2 openmm/7.7.0 ost/2.4 pdb-tools/2.5.0 promod/3.3 python/3.11.1 rmsd/1.5.1 stage/1.0.0 tmbed/1.0.0)
for mod in ${modules[@]}; do
if $(module list | grep -q "${mod}")
then
module unload ${mod}
fi
done

# required external software
module load python/3.11.1
module load biopython/1.81
module load tmbed/1.0.0
module load boost_own/1.82.0
module load openmm/7.7.0
module load ost/2.4
module load promod/3.3
module load gromacs/2023.2 gromacs=gmx
module load rmsd/1.5.1
module load pdb-tools/2.5.0

export OPENMM_CUDA_COMPILER=cuda/12.0/bin/nvcc

PATH_TO_SLURM=

INPUT_PATH=
PATH_TO_SCRIPTS=scriptsPS

PDBid=$(pwd | rev | cut -d "/" -f 1 | rev)

echo "Preparing ${PDBid} ..."
file=${INPUT_PATH}/${PDBid}/${PDBid}_protein.pdb

## Sometimes, entire chains aren't modelled in the PDB file.
## Let's assume that these chains aren't important for ligand binding in this case (If a chain isn't modelled at all, we have no checks because we have no idea where it is.)
## Therefore, let's generate FASTA files from the PDB files first and compare them to the official FASTA files in the PDB later.
# Create FASTA file
cp ${file} ${PDBid}.pdb
python3 ${PATH_TO_SCRIPTS}/pdb2fasta.py ${PDBid}.pdb

# delete folder if non-standard residues are found
if [ -f nonStandardResidueFound ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, non-standard amino acids were found in the structure. Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

## We only want to remove a protein from the data set if the ligand-binding core is a membrane protein
## We don't want to remove it due to membrane anchors etc. irrelevant to the protein-ligand interaction
## Therefore, we use our self-made FASTA file for this check
# Check whether it is a membrane protein
tmbed predict -f ${PDBid}.fasta -p ${PDBid}.pred --out-format 0 --no-use-gpu #--no-cpu-fallback

prediction=$(tail -1 ${PDBid}.pred)
if echo ${prediction} | grep -q "B" || echo ${prediction} | grep -q "b" || echo ${prediction} | grep -q "H" || echo ${prediction} | grep -q "h" || echo ${prediction} | grep -q "S"
then
echo "For ${PDBid}, some residues should be located in a membrane! Deleting folder as this protein cannot be handled by the current workflow!"
cd ..
rm -r ${PDBid}
exit 0
fi

## We have a data set of proteins co-crystallised with ligands.
## Therefore, N- or C-terminal extensions folding back and binding to the ligand should be visible in the structure
## The most probable scenario is that these extensions are either disordered or membrane anchors/tethers/etc., which we can't handle with our current workflow anyways
## Consequently, we only need to model missing residues inside the ligand-binding core unit

# Download original PDB file if our input file doesn't contain the missing residues in the sequence info
PDBid2=$(echo "$PDBid" | cut -d "_" -f 2)
set +e
wget https://www.rcsb.org/fasta/entry/${PDBid2}
set -e

if ! [ -f ${PDBid2} ]
then
echo "For ${PDBid}, no official FASTA file and thus no template sequence for gap repair could be found! Deleting folder as this protein cannot be handled by the current workflow!"
cd ..
rm -r ${PDBid}
exit 0
fi

# Collect information about gaps in the structure
pdb_delhetatm ${PDBid}.pdb > pdb_nohetatm.pdb
pdb_gap pdb_nohetatm.pdb > gaps.txt


# add missing residues to sequence of core unit
python3 ${PATH_TO_SCRIPTS}/extractSequence.py ${PDBid}.fasta ${PDBid2} ${PDBid}.pdb pdb_nohetatm.pdb gaps.txt
rm gaps.txt
rm pdb_nohetatm.pdb

# delete folder if the official PDB FASTA sequence contains the letter "X"
if [ -f xInFastaSequence ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, the official PDB FASTA sequence does not specify the identity of some residues in the sequence (X as one letter code). Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

# delete folder if the sequence of the PDB structure cannot be properly aligned to the full sequence
if [ -f alignmentError ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, the check whether residues are missing could not be executed correctly. Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

mv new.fasta ${PDBid}.fasta
rm ${PDBid2}

## We want to stay as close to the experimental structure as possible, i.e. do not remove crystal waters and ions (ligands are in a different file)
## We still have to run the protocol all the time because many structures miss just a few atoms in some side chains
## It seems that promod3 introduces very little changes to the experimental structure with the current settings (even though it rebuilds the side chains)
## We assume that the changes are so small that even side chain atoms complexated by metal ions do not move enough for the complexation to become impossible
## Strategy: accept current promod3 script, copy promod3 result back to PDB file (it's very unlikely to have ions or crystal waters resolved next to a missing loop), check RMSD to experimental structure
# fix protein structure
echo "Starting promod3"
python3 ${PATH_TO_SCRIPTS}/structureNormalisation.py ${PDBid}.pdb ${PDBid}.fasta ${PDBid}_normalised.pdb

# delete folder if the sequence of the PDB structure cannot be properly aligned to the full sequence
if [ -f normalisationError ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, promod3 didn't work correctly (alignment error). Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

# delete folder if input PDB file can't be read, e.g., because some atoms have improper naming in the original PDB
if [ -f normalisationErrorReadingPDB ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, promod3 didn't work correctly (error reading input PDB file). Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

# get missing residues from promod output (more reliable than tmbed)
tailNumber=$(grep -n "Starting promod3" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
ringPunchNumber=${tailNumber}
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "ring punch(es)")
then
ringPunchNumber=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -n "ring punch(es)" | cut -d ":" -f 1)
fi
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep -q "Resolved")
then
gaps=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep -c "Resolved")
fi
touch missingResidues.txt

for i in $(seq 1 ${gaps}); do
modelled=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Resolved" | tail -${i} | head -1)
chain=$(echo ${modelled} | cut -d "." -f 1)
chain=${chain: -1}
lowerBound=$(echo ${modelled} | cut -d "." -f 2 | cut -d "-" -f 1 | tr -d -c 0-9) # promod3 uses consecutive numbering without letters
upperBound=$(echo ${modelled} | cut -d "." -f 3 | cut -d " " -f 1 | tr -d -c 0-9) # promod3 uses consecutive numbering without letters
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep -q "Merged gap")
then
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep -q " ${chain}")
then
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep " ${chain}" | grep -q ${lowerBound})
then
lowerBound="${lowerBound} "$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep " ${chain}" | grep ${lowerBound} | cut -d "." -f 3 | cut -d " " -f 1 | tr -d -c 0-9)
upperBound=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep " ${chain}" | grep ${upperBound} | cut -d "." -f 3 | cut -d " " -f 1 | tr -d -c 0-9)" ${upperBound}"
fi
fi
fi
for h in $(seq 1 $(echo ${lowerBound} | wc -w)); do
for j in $(seq $(( $(echo ${lowerBound} | cut -d " " -f ${h}) + 1 )) $(( $(echo ${upperBound} | cut -d " " -f ${h}) - 1 ))); do
if (( ${j} == ($(echo ${upperBound} | cut -d " " -f ${h}) - 1) )) && (( ${i} == ${gaps} )) && (( ${h} == $(echo ${lowerBound} | wc -w) ))
then
echo -n "(group Protein and chain ${chain} and resnr ${j})" >> missingResidues.txt
else
echo -n "(group Protein and chain ${chain} and resnr ${j}) or " >> missingResidues.txt
fi
done
done
done

# get residues that promod3 refuses to remodel
touch residuesNotRemodelled.txt
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "incomplete backbone")
then
numberNotRemodelled=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "incomplete backbone")
for i in $(seq 1 ${numberNotRemodelled}); do
notRemodelled=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep "incomplete backbone" | head -${i} | tail -1)
currentResidueNumber=$(echo ${notRemodelled} | cut -d "." -f 2 | cut -d " " -f 1)
currentResidueNumber=${currentResidueNumber:3}
if (( ${i} == ${numberNotRemodelled} ))
then
echo -n "(group Protein and chain "$(echo ${notRemodelled} | cut -d "." -f 1)" and resname "$(echo ${notRemodelled} | cut -d "." -f 2 | cut -d " " -f 1 | cut -c1-3)" and resnr "${currentResidueNumber}")" >> residuesNotRemodelled.txt
else
echo -n "(group Protein and chain "$(echo ${notRemodelled} | cut -d "." -f 1)" and resname "$(echo ${notRemodelled} | cut -d "." -f 2 | cut -d " " -f 1 | cut -c1-3)" and resnr "${currentResidueNumber}") or "  >> residuesNotRemodelled.txt
fi
done
fi

# fuse promod3 result and crystal waters and ions, create PDB files to check RMSD
python3 ${PATH_TO_SCRIPTS}/reorganisePDBs.py ${PDBid}.pdb ${PDBid}_normalised.pdb missingResidues.txt residuesNotRemodelled.txt

# delete folder if reorganisePDBs.py fails with an index error
if [ -f indexError ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, reorganisePDBs.py failed with an index error. The root cause is most likely a subtle error in the gap identification during structure repair. Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

# delete folder if modelling is required but not successful
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Stereo-chemical problem in sidechain") && ([ -s missingResidues.txt ] || [ -s missingAtoms.txt ])
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, the crystal structure lacks atoms and could not be modelled successfully. Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

# delete folder if a mismatch between the sequence and the structure is suspected
if [ -f atomsInConectStatementsNotFound ]
then
cd ..
rm -r ${PDBid}
echo "For ${PDBid}, some atoms in CONECT statements could not be found again after modelling (ScriptError). Deleting folder as this protein should not be handled by the current workflow."
exit 0
fi

# At this stage, we know whether any residues had to be modelled or atoms in residues were missing
# If both isn't the case, go back to the original structure
# If incomplete residues had to be removed, however, take the modelled structure
if ! [ -s missingResidues.txt ]
then
if ! [ -s missingAtoms.txt ]
then
if ! [ -s residuesNotRemodelled.txt ]
then
echo "For ${PDBid}, nothing needed to be modelled. Using the original PDB file!"
# final clean up
rm ${PDBid}.fasta ${PDBid}.pred ${PDBid}_normalised.pdb original.pdb remodelled.pdb fused.pdb
mv CONECT_old.txt CONECT.txt
exit 0
fi
fi
fi

# ensure homogeneous numbering in PDB files for RMSD calculation (ignoring gaps)
gmx editconf -f original.pdb -o original2.pdb -resnr 1
mv original2.pdb original.pdb
gmx editconf -f remodelled.pdb -o remodelled2.pdb -resnr 1
mv remodelled2.pdb remodelled.pdb

# check RMSD (in units of Angstrom)
rmsd=$(calculate_rmsd original.pdb remodelled.pdb)
echo "For ${PDBid}, the RMSD between the original crystal structure and the remodelled part (in A) amounts to: ${rmsd}"
if (( $(echo "${rmsd} > 2" | bc -l) ))
then
echo "For ${PDBid}, modelling resulted in a significant deviation from the original crystal structure (more than 2 A). Deleting folder as this protein cannot be handled by the current workflow!"
cd ..
rm -r ${PDBid}
exit 0
fi
mv fused.pdb ${PDBid}.pdb
rm original.pdb remodelled.pdb

# final clean up
rm ${PDBid}.fasta ${PDBid}.pred ${PDBid}_normalised.pdb CONECT_old.txt

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload biopython/1.81
module unload tmbed/1.0.0
module unload boost_own/1.82.0
module unload openmm/7.7.0
module unload ost/2.4
module unload promod/3.3
module unload gromacs/2023.2
module unload rmsd/1.5.1
module unload pdb-tools/2.5.0

export OPENMM_CUDA_COMPILER=""
