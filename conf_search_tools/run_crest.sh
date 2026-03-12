#!/bin/bash

for i in *inp # submit crest submissions ONLY for .xyz files which have specified constraints in a .inp file.
  do
    FILE=$(basename $i .inp)
    cp crest_submission.slurm ${FILE}_crest.slurm
    sbatch ${FILE}_crest.slurm $FILE 0 ${UUFSCELL::-6}
    echo $FILE >> conformers.txt
  done
