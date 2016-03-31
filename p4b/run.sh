## solve until t=20

EXEC="./$1"

# OPT="--use-laspack NX=100 WT=13 mesh=line exact=true exo=true"
# $EXEC $OPT NT=20 dt=1 2> out/exact_err20
# $EXEC $OPT NT=40 dt=0.5 2> out/exact_err40
# $EXEC $OPT NT=80 dt=0.25 2> out/exact_err80
# $EXEC $OPT NT=160 dt=.125 2> out/exact_err160
# $EXEC $OPT NT=320 dt=0.0625 2> out/exact_err320
# $EXEC $OPT NT=640 dt=0.03125 2> out/exact_err640
# $EXEC $OPT NT=1280 dt=0.015625 2> out/exact_err1280
# $EXEC $OPT NT=2560 dt=0.0078125 2> out/exact_err2560
# $EXEC $OPT NT=5120 dt=0.0039062 2> out/exact_err5120

# OPT="--use-laspack NX=100 WT=13 mesh=line exact=false "
# $EXEC $OPT NT=20 dt=1 2> out/line_err20
# $EXEC $OPT NT=40 dt=0.5 2> out/line_err40
# $EXEC $OPT NT=80 dt=0.25 2> out/line_err80
# $EXEC $OPT NT=160 dt=.125 2> out/line_err160
# $EXEC $OPT NT=320 dt=0.0625 2> out/line_err320
# $EXEC $OPT NT=640 dt=0.03125 2> out/line_err640
# $EXEC $OPT NT=1280 dt=0.015625 2> out/line_err1280
# $EXEC $OPT NT=2560 dt=0.0078125 2> out/line_err2560
# $EXEC $OPT NT=5120 dt=0.0039062 2> out/line_err5120

OPT="--use-laspack WT=13 mesh=mesh/el.msh exact=false exo=true"
$EXEC $OPT NT=20 dt=1 2> out/el_err20
$EXEC $OPT NT=40 dt=0.5 2> out/el_err40
$EXEC $OPT NT=80 dt=0.25 2> out/el_err80
$EXEC $OPT NT=160 dt=.125 2> out/el_err160
$EXEC $OPT NT=320 dt=0.0625 2> out/el_err320
