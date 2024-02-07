#!/bin/bash
set -e

# modules consistent with GROMACS
module load LUMI/22.08
module unload cray-mpich/8.1.27 # due to tMPI
module unload cray-libsci/22.08.1.1 # due to tMPI
module load buildtools/22.08 # loads CMake/3.24.0
module load cray-python/3.9.12.1

# additional Python packages
module load typer/0.9.0

# HQ itself
module load hq/0.16.0

# start server
hq --server-dir=$(pwd)/hq-servers server start
