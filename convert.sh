#!/bin/bash

bin2vts=${HOME}/projects/bin2vts/build/bin2vts
dir=${SCRATCH}/hom_z85

for i in `ls ${dir}`; do
  ${bin2vts} -i ${dir}/${i} -dim 3 -nx 60 -ny 220 -nz 85 -o ${dir}/${i}.vts
done
