#!/bin/bash
set -e

module use /appl/local/csc/modulefiles # for access to pre-installed GROMACS versions

module load gromacs/2023.2 # for consistency

module load cray-python/3.9.12.1

# additional Python packages
module load typer/0.9.0-cpu

# HQ itself
module load hq/0.17.0-cpu

# start server
hq --server-dir=$(pwd)/hq-servers server start
