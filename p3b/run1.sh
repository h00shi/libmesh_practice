#!/bin/bash
rm out/false_dirich.txt

OPT="l=false  exact_bdry=false out=false_dirich.txt"
./readmesh_dirich wedge1 ${OPT}
./readmesh_dirich wedge2 ${OPT}
./readmesh_dirich wedge3 ${OPT}
