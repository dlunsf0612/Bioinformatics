#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o dna.2.tree.log
#SBATCH --account=dclunsford4385
#SBATCH --partition=silver

export PROJ_DIR=/export/home/bio_class/dclunsford4385/Bioinformatics
cd $PROJ_DIR

/export/share/software/biological/raxml-ng_v2.0.0/raxml-ng \
./raxml-ng --all --msa msa.fasta --model GTR --tree pars{10} --bs-trees 100

