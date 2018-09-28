#!/bin/bash
#SBATCH --image=docker:slosar/desc_lss:v0.1ut
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH -C haswell
#SBATCH --volume="/global/project/projectdirs/boss/lya/anze/work/LSST/HyperSupremeStructure-HSC-LSS/:/io"
srun -n 1 shifter /io/batch/sample_nmt.sh


