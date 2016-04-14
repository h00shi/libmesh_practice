PETSC_OPTIONS=-ksp_monitor -ksp_type cg -sub_pc_factor_mat_ordering_type qmd -sub_pc_factor_levels 2 
ROOT=../..
MPI=${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -np 
EXEC=${ROOT}/src/driver_manufactured-opt

mms1: 
	${MPI} 2 ${EXEC} mesh_name=${ROOT}/mesh/simple_box_5.msh  out_name=box_n5_o2 test_case=mms \
order=2 fe_family=1 ${PETSC_OPTIONS}
	${MPI} 2 ${EXEC} mesh_name=${ROOT}/mesh/simple_box_5.msh  out_name=box_n5_o3 test_case=mms \
order=3 fe_family=1 ${PETSC_OPTIONS}

mms2: 
	${MPI} 4 ${EXEC} mesh_name=${ROOT}/mesh/simple_box_10.msh  out_name=box_n10_o2 test_case=mms \
order=2 fe_family=1 ${PETSC_OPTIONS}
	${MPI} 4 ${EXEC} mesh_name=${ROOT}/mesh/simple_box_10.msh  out_name=box_n10_o3 test_case=mms \
order=3 fe_family=1 ${PETSC_OPTIONS}


mms3: 
	${MPI} 7 ${EXEC} mesh_name=${ROOT}/mesh/simple_box_20.msh  out_name=box_n20_o2 test_case=mms \
order=2 fe_family=1 ${PETSC_OPTIONS}
	${MPI} 7 ${EXEC} mesh_name=${ROOT}/mesh/simple_box_20.msh  out_name=box_n20_o3 test_case=mms \
order=3 fe_family=1 ${PETSC_OPTIONS}

clean:
	rm -f box*
