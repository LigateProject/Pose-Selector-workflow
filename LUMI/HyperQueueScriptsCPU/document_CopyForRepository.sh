#!/bin/bash
set -e

NUMBER=
X=

# first part
for i in $(tail -${NUMBER} HyperQueueRuns/currentlyRunning.txt | cut -d " " -f 1); do
  rm -r ${i}/rep_1/cpt ${i}/hq-work || continue
done
bash HyperQueueScriptsCPU/analyseCompletionAndCleanUp.sh >> round${X}.txt
bash HyperQueueScriptsCPU/analyseCompletionAndCleanUp.sh >> round${X}.txt

## clean up round${X}.txt and then start second part
#bash HyperQueueScriptsCPU/countPoses.sh >> round${X}.txt
#bash HyperQueueScriptsCPU/identifyRemainingPoses.sh >> round${X}.txt
#file=round${X}.txt
#number=$(grep -n "MD simulations are complete for these complexes:" ${file} | cut -d ":" -f 1)
#number=$(( ${number}+1 ))
#PDBids=$(tail -n +${number} ${file} | head -1)
#cd HyperQueueRuns
#mv startHQ startHQ_round${X}
#echo "Compressing startHQ_round${X}"
#tar czf startHQ_round${X}.tar.gz startHQ_round${X}
#rm -r startHQ_round${X} <<-eof
#y
#eof
#echo "Done"
#for PDBid in ${PDBids}; do
#  if [ -d ${PDBid} ]
#  then
#    echo "Compressing ${PDBid}"
#    tar czf ${PDBid}.tar.gz ${PDBid}
#    rm -r ${PDBid}
#    echo "Done"
#  fi
#done
#rm currentlyRunning.txt
#touch currentlyRunning.txt
#cd ..
