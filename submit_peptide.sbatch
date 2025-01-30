#!/usr/bin/env bash
#SBATCH -A <AccountInfo>
#SBATCH -p <GPUINFO>
#SBATCH --gpus=1
#SBATCH -t 00-10:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --output=./monomer-%J.out
#
# Three arguments 
# $1 - input file path
# $2 - output dir, and postfix-filename argument
# $3 - True/False for if the peptide is cyclic or not

# Define out tracking file
OUT="md_traj.csv"

# Source bashrc (as this should enable your conda env)
source $HOME/.bashrc

# conda activate
conda activate openmm

# Submit the script 
if [[ $3 == "True" ]]; then
  python md_monomer.py --file $1 --postfix-filename $2 --production --cyclic --output-dir $2 --gpu
else
  python md_monomer.py --file $1 --postfix-filename $2 --production --output-dir $2 --gpu
fi

# Write information (This isn't necessary, but might be good to do for tracking purposes)
if [[ -f $OUT ]]; then
  echo "$1,$2,`date "+%F-%T"`" >> md_traj.csv
else
  echo "filename,postfix-filename,date" > md_traj.csv
  echo "$1,$2,`date "+%F-%T"`" >> md_traj.csv
fi
