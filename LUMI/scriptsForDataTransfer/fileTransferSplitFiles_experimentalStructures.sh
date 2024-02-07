#!/bin/bash
set -e

# Data paths
## Local cluster
ROOT_PATH=/${ligate}/PSLargeScale/experimentalStructures
PATH_FILES_TO_TRANSFER=${ROOT_PATH}/inputFiles
PATH_TRANSFER_METADATA=${ROOT_PATH}/transferredSystems
PATH_TRANSFER_METADATA_DOCKING_POSES=${ROOT_PATH}/../transferredSystems

# Get a list of all systems for which experimental structures were prepared for simulation
SYSTEMS_TO_TRANSFER=($(ls ${PATH_FILES_TO_TRANSFER}))

# Info about documentation files

## ${PATH_TRANSFER_METADATA}/validationSet.txt and ${PATH_TRANSFER_METADATA_DOCKING_POSES}/validationSet.txt are not identical
## Systems listed in ${PATH_TRANSFER_METADATA_DOCKING_POSES}/validationSet.txt are all included in ${PATH_TRANSFER_METADATA_DOCKING_POSES}/toBeTransferredToMeluXina.txt
## There are seven systems in the validation set targeted for MeluXina (06571_3IQH, 10182_4HXJ, 07591_3Q8H, 06572_3IQI, 06570_3IQG, 08102_3T4P, 09830_4FGT) for which docking poses could be converted into simulation input files while the experimental structure could not (found in ${PATH_TRANSFER_METADATA_DOCKING_POSES}/validationSet.txt, but not in ${PATH_TRANSFER_METADATA}/validationSet.txt)
## There are nine systems in the validation set (01966_1SR7, 02332_1WDQ, 02495_1YC4, 02722_2ARM, 03372_2HU6, 03437_2IGV, 03438_2IGW, 04207_2R1Y, 04209_2R2B) for which the docking poses were simulated on LUMI (not found in ${PATH_TRANSFER_METADATA_DOCKING_POSES}/validationSet.txt, but in ${PATH_TRANSFER_METADATA}/validationSet.txt and ${PATH_TRANSFER_METADATA_DOCKING_POSES}/toBeTransferredToLumi.txt)

## A list of all additional experimental structures and all systems for which docking poses but not the experimental structure could be converted is available in differenceExperimentalStructuresDockingPoses (the respective lists do not distinguish between systems targeted for LUMI and MeluXina, it is one respective list for both clusters)

## copyExpStructuresToMeluXina.sh used the list of systems stored in ${PATH_TRANSFER_METADATA}/validationSet.txt, i.e. it includes the nine systems for which the docking poses were simulated on LUMI

## Given the above background info, the following selection strategy is applied:
### 1) If a system is listed in ${PATH_TRANSFER_METADATA}/validationSet.txt or ${PATH_TRANSFER_METADATA_DOCKING_POSES}/toBeTransferredToMeluXina.txt, it will be added to the list for MeluXina
### 2) Else, it will be added to the list for LUMI
### In both cases, we add a clause to make sure that systems that are already on the list are not added a second time
### => All additional experimental structures will be simulated on LUMI
### => For nine systems, docking poses simulated on LUMI will be paired with experimental structures simulated on MeluXina

# expected numbers:
## LUMI: 2601 systems with docking poses - 84 experimental structures that could not be converted - 9 experimental structures simulated on MeluXina + 13 additional experimental structures = 2521 systems
## MeluXina: 3128 (total) - 2521 (LUMI) = 607 systems = 610 (currently prepared docking poses) + 9 (experimental structures for LUMI systems) - 12 (experimental structures that could not be converted)

for PDBid in ${SYSTEMS_TO_TRANSFER[@]}; do
  if (grep -q ${PDBid} ${PATH_TRANSFER_METADATA}/validationSet.txt || grep -q ${PDBid} ${PATH_TRANSFER_METADATA_DOCKING_POSES}/toBeTransferredToMeluXina.txt) && (! grep -q ${PDBid} ${PATH_TRANSFER_METADATA}/toBeTransferredToMeluXina.txt)
  then
    echo ${PDBid} >> ${PATH_TRANSFER_METADATA}/toBeTransferredToMeluXina.txt
  elif ! grep -q ${PDBid} ${PATH_TRANSFER_METADATA}/toBeTransferredToLumi.txt
  then
    echo ${PDBid} >> ${PATH_TRANSFER_METADATA}/toBeTransferredToLumi.txt
  fi
done
