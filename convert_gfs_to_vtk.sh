#!/bin/bash

for i in `seq 0 1 0`; do
for j in `seq 0 1 9`; do
for k in `seq 0 1 9`; do

sleep 1
echo 't='${i}${j}${k}
gerris2D -e "OutputSimulation {istep = 1} snap-0.${i}${j}${k}.vtk {format=VTK}" snap-0.${i}${j}${k}.gfs > a

done
done
done
done
done
