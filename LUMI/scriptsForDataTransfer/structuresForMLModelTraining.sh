#!/bin/bash

#SBATCH --job-name=PreparePSDataForLongTermStorage
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH -t 24:00:00
#SBATCH -o slurm-%J.out

set -e

# load required modules
module load forSTaGE/openbabel/3.1.1
module load gromacs/2023.2

# Paths
cwd=/scratch
PATH_ARCHIVES=${ligate}/PSLargeScale/archivedSystems
PATH_TO_SCRIPTS=LUMI/PostProcessing

# Data
cd ${PATH_ARCHIVES}
ARCHIVES=($(ls))

# Select batch of archives to work on
BATCH_SIZE=10
START_INDEX=2591 # start index for this batch

# Go to the node's local hard drive not to overload the file system
mkdir -p ${cwd}
cd ${cwd}

# Work on archives
index=${START_INDEX}

while (( ${index} < $(( ${BATCH_SIZE}+${START_INDEX} )) )); do

  if (( ${index}=${#ARCHIVES[@]} ))
  then
    break
  fi

  folder=$(echo ${ARCHIVES[${index}]} | cut -d "." -f 1)
  echo "Starting to generate ML model input structure files for ${folder}!"

  # copy relevant data
  cp ${PATH_ARCHIVES}/${ARCHIVES[${index}]} .

  # we should only produce ML model input files for poses that could be simulated successfully
  if [ -d ${ps1}/${folder} ]
  then
    directory=${ps1}
  elif [ -d ${ps2}/${folder} ]
  then
    directory=${ps2}
  elif [ -d ${ps3}/${folder} ]
  then
    directory=${ps3}
  fi

  # inflate
  tar xzf ${ARCHIVES[${index}]}
  rm ${ARCHIVES[${index}]}
  mv ${folder} ${folder}_archive

  # produce input structure file for ML model
  mkdir ${folder}
  cd ${folder}
  for pose in $(ls ${directory}/${folder} | grep "pose"); do
    pose_folder=$(echo ${pose} | cut -d "." -f 1)
    mkdir ${pose_folder}

    # repair ions.gro if the structure is broken across periodic boundaries
    lineNumber=$(grep -n "SOL" ../${folder}_archive/${pose_folder}/ions.gro | head -1 | cut -d ":" -f 1)
    head -$(( ${lineNumber}-1 )) ../${folder}_archive/${pose_folder}/ions.gro >> input.gro
    tail -1 ../${folder}_archive/${pose_folder}/ions.gro >> input.gro
    python3 ${PATH_TO_SCRIPTS}/makeStartingStructureWhole.py
    ## remove last line of output.gro
    sed -i '$ d' output.gro
    tail -n +${lineNumber} ../${folder}_archive/${pose_folder}/ions.gro >> output.gro
    mv output.gro ${pose_folder}/ions.gro
    rm input.gro intermediary.gro

    if ! [ -f indexRerun.ndx ]
    then 
      # identify the correct keyword to grep for a potential peptide ligand
      ligandIdentifier=$(grep "Protein" ../${folder}_archive/topol_amber_${pose_folder}.top | tail -1 | cut -d " " -f 1)

      # create new index file
      ## system may be charged, and there might be a harmless mismatch in the name of hydrogen atoms
      gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c ${pose_folder}/ions.gro -p ../${folder}_archive/topol_amber_${pose_folder}.top -o selection.tpr -po selectionOut.mdp -maxwarn 2
      ## identify ligand correctly
      ### small organic molecule or any combinations or special cases
      val1=$(grep -n "; Compound" ../${folder}_archive/topol_amber_${pose_folder}.top | cut -d ":" -f 1)
      if grep -q "MOL " ../${folder}_archive/topol_amber_${pose_folder}.top
      then
        val2=$(grep -n "MOL " ../${folder}_archive/topol_amber_${pose_folder}.top | cut -d ":" -f 1)
      ### peptide ligand
      elif grep -q "${ligandIdentifier} " ../${folder}_archive/topol_amber_${pose_folder}.top
        then
        val2=$(grep -n "${ligandIdentifier} " ../${folder}_archive/topol_amber_${pose_folder}.top | cut -d ":" -f 1)
      ### TODO add clauses for DNA and RNA ligands
      fi
      diff=$(( ${val2} - ${val1} ))
      
      ## create separate bash script to run gmx select
      ligand='"ligand" molecule '${diff}''
      ligand=\'${ligand}\'
      protein='"protein" group Protein and not molecule '${diff}''
      protein=\'${protein}\'
      
      gmx_command="gmx select -s selection.tpr -on indexRerun.ndx -select ${ligand} ${protein}"
      
      echo "#!/bin/bash" >> createIndex.sh
      echo "set -e" >> createIndex.sh
      echo "" >> createIndex.sh
      echo $gmx_command >> createIndex.sh
      
      ## actually create index file
      bash createIndex.sh
      
      ## clean up
      rm createIndex.sh selection.tpr selectionOut.mdp
    fi

    # create required structure files with GROMACS
    cd ${pose_folder}
    gmx trjconv -f ions.gro -s ../../${folder}_archive/${pose_folder}/simulation.tpr -n ../indexRerun.ndx -o protein.pdb <<-eof
protein
eof
    gmx trjconv -f ions.gro -s ../../${folder}_archive/${pose_folder}/simulation.tpr -n ../indexRerun.ndx -o ligand.gro <<-eof
ligand
eof

    # convert ligand file to MOL2
    obabel "ligand.gro" -O "ligand.mol2"

    # clean up
    rm ligand.gro ions.gro
    cd ..

  done

  if [ -f indexRerun.ndx ]
  then
    # clean up
    rm indexRerun.ndx
    cd ..

    # compress
    tar czf ${ARCHIVES[${index}]} ${folder}

    # copy archive to final destination
    cp ${ARCHIVES[${index}]} $pslt/MLStructureFiles

  else
    cd ..
  fi

  # clean up
  if [ -f ${ARCHIVES[${index}]} ]
  then
    rm ${ARCHIVES[${index}]}
  fi
  rm -r ${folder} ${folder}_archive

  index=$(( ${index}+1 ))
  echo "Done generating ML model input structure files for ${folder}!"

done

# Completed
echo "Complete!"
