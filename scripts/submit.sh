#!/bin/sh


ROOT="/tcghome/ygeorgi/fock"
LD_LIBRARY_PATH="$ROOT/gcc/lib64:$ROOT/lib64"
PATH="/bin"
export LD_LIBRARY_PATH PATH

MPIRUN="$ROOT/mpi/bin/mpirun"
EXE="$ROOT/rotd/multi"

MACHINEFILE=""
OPTIONS=""

OMP_NUM_THREADS="1"
 
while test "x$1" != x ; do
 case $1 in
  -f | --machinefile 			) shift; MACHINEFILE="$1"	;;
  -t | --threads			) shift; OMP_NUM_THREADS="$1"	;;
  *					) OPTIONS="$OPTIONS $1" 	;;
 esac
 shift
done

if [ "x$MACHINEFILE" == "x" -o "x$OPTIONS" == "x" ]; then
  echo "usage: `basename $0` {-f|--machinefile} <nodes_file_name> {-t | --threads} <threads_no> input_file_name"
  exit 0
fi

if [ "x$OMP_NUM_THREADS" == "x" ]; then
	exec $MPIRUN -x LD_LIBRARY_PATH -machinefile $MACHINEFILE $EXE $OPTIONS
else
	export OMP_NUM_THREADS
	exec $MPIRUN -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -machinefile $MACHINEFILE $EXE $OPTIONS
fi
