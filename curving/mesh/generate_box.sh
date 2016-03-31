#!/bin/bash


ID="0 1 2 3"
NLAYER=(5 10 20 40)

for ii in $ID
do
	echo "Generating simple_box for NLAYER=${NLAYER[ii]}"
	gmsh simple_box.geo -3 -o "simple_box_${NLAYER[ii]}.msh" -setnumber n_layer ${NLAYER[ii]}
done

