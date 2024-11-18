#!/bin/bash
#SBATCH --job-name=PEG_conc
#SBATCH -t 4-04:00:00 
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p xeon-p8
#SBATCH --mem=16GB

g++ -std=c++11 Module_4_5/Module_4_5.cpp -o nerd
./nerd 0.018 40000
