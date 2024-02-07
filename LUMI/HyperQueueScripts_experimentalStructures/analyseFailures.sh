#!/bin/bash
set -e

NUMBER=

# identify failed runs by the presence of a cpt folder 
find -name cpt >> failed.txt
allFailures=$(cat failed.txt | wc -l)
echo "Identified all potential failures"

# exclude failures due to dumped CPU cores
toBeAnalysedFurther=$(for i in $(seq 1 ${allFailures}); do crashed=$(grep "core dumped" $(head -${i} failed.txt | tail -1)/../../hq-work/job-*/stderr) || echo ${i}; done)
echo "$(( ${allFailures}-$(echo ${toBeAnalysedFurther} | wc -w) )) CPU core dumps"
allFailures=$(echo ${toBeAnalysedFurther} | wc -w)

# exclude failures due to additional segmentation faults
toBeAnalysedFurther=$(for i in ${toBeAnalysedFurther}; do memoryLeak=$(grep "Segmentation fault" $(head -${i} failed.txt | tail -1)/../../hq-work/job-*/stderr) || echo ${i}; done)
echo "$(( ${allFailures}-$(echo ${toBeAnalysedFurther} | wc -w) )) additional segmentation faults"
allFailures=$(echo ${toBeAnalysedFurther} | wc -w)

# exclude failures due to GPU memory access problems
toBeAnalysedFurther=$(for i in ${toBeAnalysedFurther}; do memoryLeak=$(grep "Reason: Unknown" $(head -${i} failed.txt | tail -1)/../../hq-work/job-*/stderr) || echo ${i}; done)
echo "$(( ${allFailures}-$(echo ${toBeAnalysedFurther} | wc -w) )) additional GPU memory leaks"
allFailures=$(echo ${toBeAnalysedFurther} | wc -w)

# exclude failures due to mkdir errors
backup=${toBeAnalysedFurther}
toBeAnalysedFurther=$(for i in ${toBeAnalysedFurther}; do mkdirErrorsDummy=$(grep "mkdir: cannot" $(head -${i} failed.txt | tail -1)/../../hq-work/job-*/stderr) || echo ${i}; done)
echo "$(( ${allFailures}-$(echo ${toBeAnalysedFurther} | wc -w) )) errors due to mkdir for existing directory"
allFailures=$(echo ${toBeAnalysedFurther} | wc -w)
# identify all mkdir errors
mkdirErrors=$(for i in ${backup}; do dummy=$(echo ${toBeAnalysedFurther} | grep -w "${i}") || echo ${i}; done)
for i in ${mkdirErrors}; do
  rm -r $(head -${i} failed.txt | tail -1)
done
echo "Removed the corresponding $(echo ${mkdirErrors} | wc -w) cpt folders."

# ignore poses that will be dealt with by document.sh
toBeAnalysedFurther=$(for i in ${toBeAnalysedFurther}; do pose=$(head -${i} failed.txt | tail -1); document=$(tail -${NUMBER} HyperQueueRuns_experimentalStructures/currentlyRunning.txt | grep "${pose::-10}") || echo ${i}; done)
echo "$(( ${allFailures}-$(echo ${toBeAnalysedFurther} | wc -w) )) poses are handled by document.sh"
allFailures=$(echo ${toBeAnalysedFurther} | wc -w)

# did we cover everything?
echo "Residual failures: ${allFailures}"
echo "Have a look at:"
for i in ${toBeAnalysedFurther}; do
  echo "less $(head -${i} failed.txt | tail -1)/../../hq-work/job-*/stderr"
done
echo "Or if the simulations were just stopped prematurely, run:"
for i in ${toBeAnalysedFurther}; do
  echo "rm -r $(head -${i} failed.txt | tail -1)"
done

# clean up
rm failed.txt
