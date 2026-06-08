#!/bin/bash

INPUTFILE="protein"
OUTPUTDIR="${INPUTFILE}_sd"
RANDOMSEED=0
NUM_STRUCTURES=100

export PATH="${HOME}/localcolabfold/localcolabfold/colabfold-conda/bin:$PATH"

# Shallow MSA, 3 recycles, with dropout, 100 structures
colabfold_batch \
  --model-type alphafold2_ptm \
  --num-recycle 3 \
  --num-models 1 \
  --recycle-early-stop-tolerance 0.5 \
  --max-seq 16 \
  --max-extra-seq 32 \
  --random-seed ${RANDOMSEED} \
  --num-seeds ${NUM_STRUCTURES} \
  --num-relax ${NUM_STRUCTURES} \
  --use-dropout \
  --amber \
  --templates \
  --use-gpu-relax \
  ${INPUTFILE}.fasta \
  ${OUTPUTDIR}
