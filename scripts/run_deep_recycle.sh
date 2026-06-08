#!/bin/bash

INPUTFILE="protein"
OUTPUTDIR="${INPUTFILE}_dr"
RANDOMSEED=0
NUM_STRUCTURES=100

export PATH="${HOME}/localcolabfold/localcolabfold/colabfold-conda/bin:$PATH"

# Default MSA, 20 recycles, 100 structures
colabfold_batch \
  --model-type alphafold2_ptm \
  --num-recycle 20 \
  --num-models 1 \
  --recycle-early-stop-tolerance 0.5 \
  --random-seed ${RANDOMSEED} \
  --num-seeds ${NUM_STRUCTURES} \
  --num-relax ${NUM_STRUCTURES} \
  --amber \
  --templates \
  --use-gpu-relax \
  ${INPUTFILE}.fasta \
  ${OUTPUTDIR}
