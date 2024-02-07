#!/bin/bash -l
set -e

# Load all necessary modules and prepare your environment here
# hwloc: /usr/bin/hwloc-*; version 2.6.0

module use /appl/local/csc/modulefiles # for hipSYCL

#module load LUMI/22.08
module load hipsycl/0.9.4 # loads its own rocm
#module unload cray-mpich/8.1.27 # due to tMPI
#module unload cray-libsci/22.08.1.1 # due to tMPI
#module load buildtools/22.08 # loads CMake/3.24.0
module load cray-fftw/3.3.10.1
#module load cray-python/3.9.12.1

export LD_LIBRARY_PATH=/appl/lumi/SW/LUMI-22.08/G/EB/rocm/5.3.3/llvm/lib:$LD_LIBRARY_PATH # to be on the safe side: it's not in the module environment but needed by the clang compiler

module load gromacs/2023.2-sycl-tmpi

# set environment variables
NTMPI=1
NTOMP=6

export OMP_NUM_THREADS=${NTOMP}

if [[ -z "${NUM_ITERS}" ]]; then
  echo "NUM_ITERS has to be specified"
  exit 1
fi

# CPU-GPU mapping on LUMI
#    CPU cores       GPU id  virtual CPU cores
#   0-7,  8-15     2 (4, 5)     64-71,   72-79
# 16-23, 24-31     1 (2, 3)     80-87,   88-95
# 32-39, 40-47     3 (6, 7)    96-103, 104-111
# 48-55, 56-63     0 (0, 1)   112-119, 120-127

# different core numberings (actual core ID is the one used by the operating system)
# GROMACS pinning ID:   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27
# actual core ID:       1   65    2   66    3   67    4   68    5   69    6   70    7   71    9   73   10   74   11   75   12   76   13   77   14   78   15   79
# hwloc:                0    0    1    1    2    2    3    3    4    4    5    5    6    6    7    7    8    8    9    9   10   10   11   11   12   12   13   13

# GROMACS pinning ID:  28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54   55
# actual core ID:      17   81   18   82   19   83   20   84   21   85   22   86   23   87   25   89   26   90   27   91   28   92   29   93   30   94   31   95
# hwloc:               14   14   15   15   16   16   17   17   18   18   19   19   20   20   21   21   22   22   23   23   24   24   25   25   26   26   27   27

# GROMACS pinning ID:  56   57   58   59   60   61   62   63   64   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80   81   82   83
# actual core ID:      33   97   34   98   35   99   36  100   37  101   38  102   39  103   41  105   42  106   43  107   44  108   45  109   46  110   47  111
# hwloc:               28   28   29   29   30   30   31   31   32   32   33   33   34   34   35   35   36   36   37   37   38   38   39   39   40   40   41   41

# GROMACS pinning ID:  84   85   86   87   88   89   90   91   92   93   94   95   96   97   98   99  100  101  102  103  104  105  106  107  108  109  110  111
# actual core ID:      49  113   50  114   51  115   52  116   53  117   54  118   55  119   57  121   58  122   59  123   60  124   61  125   62  126   63  127
# hwloc:               42   42   43   43   44   44   45   45   46   46   47   47   48   48   49   49   50   50   51   51   52   52   53   53   54   54   55   55

# hwloc considers only physical cores, does not see cores blocked by the operating system and uses a consecutive numbering (not GROMACS's numbering)
HWLOCBEGIN=(dummy 0 14 28 42 7 21 35 49)
HWLOCEND=(dummy 6 20 34 48 13 27 41 55)
# make sure that GROMACS threads do not end up on the same core as SYCL worker threads
PINOFFSET=(dummy 2 30 58 86 16 44 72 100)
GPUID=(dummy 4 2 6 0 5 3 7 1)

# run simulations
pids=()

# double-check HyperQueue
echo ${CWD} $(date) $(hostname) >> HyperQueueRuns_experimentalStructures/currentlyRunning.txt

# work in the nodes memory /tmp
cd /tmp
cp -r ${CWD} .
archive=$(echo ${CWD} | rev | cut -d "/" -f 1 | rev)
cd ${archive}
# signal to analysis script that the job was started
mkdir -p ${CWD}/rep_1/cpt

for (( i=1; i<=${NUM_ITERS}; i++ )); do
    cd rep_${i}
    mkdir -p cpt

    ROCR_VISIBLE_DEVICES=${GPUID[${i}]} OMP_PLACES=cores  OMP_PROC_BIND=close hwloc-bind --cpubind "core:${HWLOCBEGIN[${i}]}-${HWLOCEND[${i}]}" \
    gmx mdrun  -pin on -pinoffset ${PINOFFSET[${i}]} -pinstride 1 \
               -deffnm simulation \
               -cpnum -cpt 60 -cpo cpt/state \
               -ntmpi ${NTMPI} -ntomp ${NTOMP} \
               -maxh 1 \
               -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu &

    # Store the PID for later
    pids+=($!)
    cd ..
done

# Wait for all pids
# This approach makes sure that wait will fail if any of the children have failed
for pid in ${pids[@]}; do
    wait $pid
done

# clean up to reduce the number of files (only for successfully run simulations)
if (( $(find rep_* -name simulation.gro | wc -l) == ${NUM_ITERS} ))
then
    # remove all GROMACS files but the .tpr and .xtc
    rm -r rep_*/cpt/ rep_*/simulation.edr rep_*/simulation.gro rep_*/simulation.log
    rm -r hq-work
    # HQ needs the output file to be located on the file system
    cp ${CWD}/hq-work/job-*/stderr .
    # convert pose into tar archive
    cd ..
    tar czf ${archive}.tar.gz ${archive}
    cp -r ${archive}.tar.gz ${CWD}/..
    rm -r ${CWD}
    rm -r ${archive} ${archive}.tar.gz
fi
