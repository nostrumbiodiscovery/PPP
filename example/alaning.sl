#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=72
#SBATCH --mem-per-cpu=1000

module purge

#2018 GNU Alaning Scanning
module load GCCcore/6.3.0 Python/2.7.14-foss-2018a
module load Python/2.7.10-foss-2018a Boost/1.66.0-foss-2018a JsonCpp/1.8.4-foss-2018a Crypto++/6.1.0-intel-2018a patchelf/0.9-foss-2018a wjelement/1.3-foss-2018a CMake/3.7.2-GCCcore-6.3.0 GCCcore/6.4.0 OpenMPI
export PYTHONPATH=/home/dsoler/:/home/dsoler/new/AdaptivePELE/:$PYTHONPATH
python ../mutate_PELE.py chC_3ngb_plopped.pdb mutations.txt --box -37.866  86.965 -49.109 --cpus 72

 
