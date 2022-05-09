#!/bin/bash
#SBATCH -J small_bert
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus-gpu
#SBATCH --qos=campus-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=0-04:00:00
#SBATCH --error=small_bert.e%J
#SBATCH --output=small_bert.o%J

cd /lustre/isaac/scratch/oqueen/CodonBert/LargeMLM/ranks
conda activate /lustre/isaac/scratch/oqueen/codonbert

python3 train_roberta_small.py