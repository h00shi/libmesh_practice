#!/bin/bash
rm out/true_dirich.txt

OPT="l=true  exact_bdry=false out=true_dirich.txt"
./readmesh_dirich wedge1 ${OPT}
./readmesh_dirich wedge2 ${OPT}
./readmesh_dirich wedge3 ${OPT}
