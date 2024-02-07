#!/bin/bash
set -e; shopt -s expand_aliases

script=runInputFileGeneration.sh

PATH_TO_SCRIPTS=$(pwd)
cd $(lsvol | cut -d " " -f 4)/PSLargeScale
mkdir inputFiles

# final value of sequence is 18441
for i in $(seq 1 10 18441); do
cp ${PATH_TO_SCRIPTS}/${script} ${script}_${i}
sed -i ${script}_${i} -e "s/#tbr"${i}"PDBids/PDBids/g"
sbatch ${script}_${i}
done
