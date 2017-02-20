#!/bin/bash
ID="0 1 2 3"
NLAYER=(10 20 30 40)

for ii in $ID
do
	echo "Generating sphere for NLAYER=${NLAYER[ii]}"
	gmsh sphere.geo -3 -o "sphere_${NLAYER[ii]}.msh" -setnumber n_r ${NLAYER[ii]}
done
