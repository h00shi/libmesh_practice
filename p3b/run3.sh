#!/bin/bash
rm out/true_ex.txt

OPT="l=true exact_bdry=true out=true_ex.txt"
./readmesh_dirich wedge1 ${OPT}
./readmesh_dirich wedge2 ${OPT}
./readmesh_dirich wedge3 ${OPT}
