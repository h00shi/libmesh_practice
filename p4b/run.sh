## solve until t=20

EXEC="./$1"
OPT="--use-laspack NX=100 WT=13"

$EXEC $OPT NT=20 dt=1 2> out/err20
$EXEC $OPT NT=40 dt=0.5 2> out/err40
$EXEC $OPT NT=80 dt=0.25 2> out/err80
$EXEC $OPT NT=160 dt=.125 2> out/err160
$EXEC $OPT NT=320 dt=0.0625 2> out/err320
$EXEC $OPT NT=640 dt=0.03125 2> out/err640
$EXEC $OPT NT=1280 dt=0.015625 2> out/err1280
$EXEC $OPT NT=2560 dt=0.0078125 2> out/err2560
$EXEC $OPT NT=5120 dt=0.0039062 2> out/err5120

