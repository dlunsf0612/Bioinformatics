#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o tree.test.log
#SBATCH --account=dclunsford4385
#SBATCH --partition=silver

module load biological/java

export PROJ_DIR=/export/home/bio_class/dclunsford4385/Bioinformatics
cd $PROJ_DIR

/export/share/software/biological/raxml-ng_v2.0.0/raxml-ng \
./raxml-ng --all --msa msa.fasta --model LG+G8+F --tree pars{10} --bs-trees 100

