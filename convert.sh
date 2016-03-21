#!/bin/bash

bin2vts=${HOME}/bin2vts/build/bin2vts

#for i in `seq 1 2`; do
for i in `ls output`; do
  #${bin2vts} -i output/${i} -dim 3 -nx 60 -ny 220 -nz 1 -o output/${i}.vts
  ${bin2vts} -i output/${i} -dim 2 -nx 60 -ny 220 -o output/${i}.vts
done
