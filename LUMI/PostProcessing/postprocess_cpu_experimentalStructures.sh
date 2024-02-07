#!/bin/bash -l
#SBATCH --job-name=post-process
#SBATCH --output=PostProcessing.o%j
#SBATCH --error=PostProcessing.e%j
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-23:59:59
set -e

module use /appl/local/csc/modulefiles # for access to pre-installed GROMACS versions

module load gromacs/2023.2

inputPath=
outputPath=
projectPath=

# we copy the TAR archives with the topology information before simulating the system
upperLimit=2521
numberOfSystems=500
PDBids=($(ls ${inputPath} | head -${upperLimit} | tail -${numberOfSystems}))
echo "Range: $(( ${upperLimit}-${numberOfSystems}+1 )) to ${upperLimit}"

PATH_TO_SCRIPTS=${projectPath}/PostProcessing
solvationModel="gb"

wd="/tmp"

for PDBid in ${PDBids[@]}; do
  # check if poses of the system have to be post-processed
  ## skip if the system is not ready for post-processing
  if ! [ -f ${outputPath}/${PDBid}/${PDBid}.tar.gz ]
  then
    continue
  ## skip if another instance is already working on the system
  elif [ -f ${outputPath}/${PDBid}/backup.tar.gz ]
  then
    continue
  ## skip if the system has been post-processed completely
  elif [ -f ${outputPath}/${PDBid}/readyForAnalysis ]
  then
    continue
  fi

  ## skip if the system has not been simulated completely
  if (( $(ls ${inputPath}/${PDBid}/pose* | wc -w)>$(ls ${inputPath}/${PDBid}/*.tar.gz | wc -w) ))
  then
    continue
  fi

  ## skip if simulations have failed for all poses of the system
  if (( $(ls ${inputPath}/${PDBid}/pose* | wc -w)==0 ))
  then
    continue
  fi

  # create backup to be able to re-start in case of failure
  cp ${outputPath}/${PDBid}/${PDBid}.tar.gz ${outputPath}/${PDBid}/backup.tar.gz

  # get relevant data
  cd ${wd}
  mkdir ${PDBid}
  cd ${PDBid}
  cp ${outputPath}/${PDBid}/${PDBid}.tar.gz .
  tar xzf ${PDBid}.tar.gz
  rm ${PDBid}.tar.gz

  # identify the correct keyword to grep for a potential peptide ligand
  ligandIdentifier=$(grep "Protein" ${PDBid}/topol_amber.top | tail -1 | cut -d " " -f 1)
  GROFile=$(ls ${PDBid}/ions* | head -1)

  # create new index file
  # can be done here because we enforce a consistent protonation pattern among poses
  ## system may be charged, and there might be a harmless mismatch in the name of hydrogen atoms
  gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c ${GROFile} -p ${PDBid}/topol_amber.top -o selection.tpr -po selectionOut.mdp -maxwarn 2
  ## identify ligand correctly
  ### small organic molecule or any combinations or special cases
  val1=$(grep -n "; Compound" ${PDBid}/topol_amber.top | cut -d ":" -f 1)
  if grep -q "MOL " ${PDBid}/topol_amber.top
  then
    val2=$(grep -n "MOL " ${PDBid}/topol_amber.top | cut -d ":" -f 1)
  ### peptide ligand
  elif grep -q "${ligandIdentifier} " ${PDBid}/topol_amber.top
    then
    val2=$(grep -n "${ligandIdentifier} " ${PDBid}/topol_amber.top | cut -d ":" -f 1)
  ### TODO add clauses for DNA and RNA ligands
  fi
  diff=$(( ${val2} - ${val1} ))

  ## create separate bash script to run gmx select
  ligand='"ligand" molecule '${diff}''
  ligand=\'${ligand}\'
  protein='"protein" not molecule '${diff}''
  protein=\'${protein}\'
  complex='"complex" group System'
  complex=\'${complex}\'
  complexWithoutIons='"complex" not group Ion'
  complexWithoutIons=\'${complexWithoutIons}\'

  gmx_command="gmx select -s selection.tpr -on indexRerun.ndx -select ${ligand} ${protein} ${complex}"

  echo "#!/bin/bash" >> createIndex.sh
  echo "set -e" >> createIndex.sh
  echo "" >> createIndex.sh
  echo $gmx_command >> createIndex.sh

  ## actually create index file
  bash createIndex.sh

  ## clean up
  rm createIndex.sh selection.tpr selectionOut.mdp

  for h in $(ls ${inputPath}/${PDBid} -I "fail*"); do
    pose=$(echo ${h} | cut -d "." -f 1)
    # get relevant data
    cp ${inputPath}/${PDBid}/${pose}.tar.gz .
    tar xzf ${pose}.tar.gz
    rm ${pose}.tar.gz

    # make sure that the structure used to generate the .tpr file is not broken across periodic boundary conditions
    cp ${PDBid}/ions_${pose}.gro input.gro
    python3 ${PATH_TO_SCRIPTS}/makeStartingStructureWhole.py
    mv output.gro ${PDBid}/ions_${pose}.gro
    rm input.gro intermediary.gro

    cd ${pose}

    ## grompp
    gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/simulationAnalysis.mdp -c ../${PDBid}/ions_${pose}.gro -p ../${PDBid}/topol_amber.top -n ../indexRerun.ndx -o rerun.tpr -po rerunOut.mdp -maxwarn 2

    for i in $(seq -w 1 8); do
      cd rep_${i}

      # prepare trajectory for binding free energy calculation
      ## remove periodic boundary conditions from trajectory
      gmx trjconv -f simulation.xtc -s ../rerun.tpr -n ../../indexRerun.ndx -o ${solvationModel}_1.xtc -b 100 -e 100 -pbc nojump <<-eof
complex
complex
eof
      gmx trjconv -f ${solvationModel}_1.xtc -s ../rerun.tpr -n ../../indexRerun.ndx -o ${solvationModel}.xtc -pbc no -fit progressive <<-eof
complex
complex
eof
      # clean up
      rm ${solvationModel}_1.xtc simulation.xtc
      cd ..
    done

    # remove ions complexated by the protein
    if grep -q "Ion" ../${PDBid}/topol_amber.top
    then
      # obtain .gro file without complexated ions
        gmx trjconv -f rep_1/${solvationModel}.xtc -s rerun.tpr -n ../indexRerun.ndx -o noIons.gro -dump 0 <<-eof
complexWithoutIons
eof
      # remove ions from trajectories
      for i in $(seq -w 1 8); do
        cd rep_${i}
        gmx trjconv -f ${solvationModel}.xtc -s ../rerun.tpr -n ../../indexRerun.ndx -o ${solvationModel}_noIons.xtc <<-eof
complexWithoutIons
eof
        # clean up
        mv ${solvationModel}_noIons.xtc ${solvationModel}.xtc
        cd ..
      done
      # update topology
      cp ../${PDBid}/topol_amber.top .
      sed -i "/Ion/d" topol_amber.top
      # create new index file
      ## system may be charged
      gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c noIons.gro -p topol_amber.top -o selection.tpr -po selectionOut.mdp -maxwarn 1
      ## identify ligand correctly
      ### small organic molecule or any combinations or special cases
      val1=$(grep -n "; Compound" topol_amber.top | cut -d ":" -f 1)
      if grep -q "MOL " topol_amber.top
      then
        val2=$(grep -n "MOL " topol_amber.top | cut -d ":" -f 1)
      ### peptide ligand
      elif grep -q "${ligandIdentifier} " topol_amber.top
        then
        val2=$(grep -n "${ligandIdentifier} " topol_amber.top | cut -d ":" -f 1)
      ### TODO add clauses for DNA and RNA ligands
      fi
      diff=$(( ${val2} - ${val1} ))

      ## create separate bash script to run gmx select
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

      ## actually create index file
      bash createIndex.sh

      ## clean up
      rm createIndex.sh selection.tpr selectionOut.mdp
      mv noIons.gro ../${PDBid}/ions_${pose}.gro
      mv topol_amber.top ../${PDBid}/topol_amber_new.top
      mv indexRerun.ndx ../indexRerun_new.ndx

      ## update TPR file
      rm rerun.tpr rerunOut.mdp
      gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/simulationAnalysis.mdp -c ../${PDBid}/ions_${pose}.gro -p ../${PDBid}/topol_amber_new.top -n ../indexRerun_new.ndx -o rerun.tpr -po rerunOut.mdp -maxwarn 1
    fi


    # clean up
    mv rerun.tpr gb.tpr
    rm rerunOut.mdp
    if [ -f stderr ]
    then
      rm stderr
    fi
    cd ..
    # make sure Akash has the necessary permissions
    chmod -R go+rwx ${pose}
    # compress
    tar czf ${pose}.tar.gz ${pose}
    # copy post-processed pose to file system
    cp ${pose}.tar.gz ${outputPath}/${PDBid}
    # clean up
    rm -r ${pose} ${pose}.tar.gz
  done

  # clean up
  rm ${PDBid}/ions_pose_*
  if [ -f indexRerun_new.ndx ]
  then
    mv indexRerun_new.ndx indexRerun.ndx
  fi
  if [ -f ${PDBid}/topol_amber_new.top ]
  then
    mv ${PDBid}/topol_amber_new.top ${PDBid}/topol_amber.top
  fi
  # add index file
  cp indexRerun.ndx ${PDBid}
  # make sure Akash has the necessary permissions
  chmod -R go+rwx ${PDBid}
  # compress and copy
  tar czf ${PDBid}.tar.gz ${PDBid}
  cp ${PDBid}.tar.gz ${outputPath}/${PDBid}
  # clean up
  cd ..
  rm -r ${PDBid}
  # make sure Akash has the permission to read and edit the tar files
  chmod o+rw ${outputPath}/${PDBid}/*.tar.gz

  # mark system as done
  touch ${outputPath}/${PDBid}/readyForAnalysis
  rm ${outputPath}/${PDBid}/backup.tar.gz
done
