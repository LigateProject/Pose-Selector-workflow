## Execution
1) Copy input files from local cluster to LUMI
2) Prepare TPR files using ../prepareFolderStructure.sh
3) Modify list of system in runner.sh and start jobs

## Sequence of steps for re-starts:
1) Check that last TPR file generation was completed successfully (grep -c "Done" PSFolderPreparation.o%j) and remove corresponding SLURM files
2) Shut down HQ properly
3) run analyseFailures.sh and document.sh
4) Re-start HyperQueue run
5) Copy over new systems
6) Start next round of TPR file generation
