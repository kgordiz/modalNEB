#!/bin/bash
#SBATCH -N 4
##SBATCH -n 40
#SBATCH --time=72:00:00
#SBATCH --constraint=centos7
#SBATCH --mem=186G
#SBATCH --ntasks-per-node=40
#SBATCH -p sched_mit_ase
#SBATCH -o slurm.out 
#SBATCH -e slurm.err
#SBATCH -J pythonw

module purge
module load gcc/6.3.0
module load openmpi/gcc/64/1.8.1
module load jdk/1.8.0_121
module load cmake/3.9.6
module load boost/1.70.0

eval "$(conda shell.bash hook)"
conda activate /home/kgordiz/miniconda2/envs/Kia1/
sleep 10
python code0_main_with_loop.py
conda deactivate
