#!/bin/bash

ID="0 1 2 3"
NLAYER=(8 16 32 64)

for ii in $ID
do
	echo "Generating bump for NLAYER=${NLAYER[ii]}"
	gmsh bump.geo -3 -o "bump_${NLAYER[ii]}.msh" -setnumber nx1 ${NLAYER[ii]}
done
