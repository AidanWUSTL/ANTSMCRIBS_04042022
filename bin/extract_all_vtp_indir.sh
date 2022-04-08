#!/bin/bash

if [ ! -f "$1" -o ! -d "$2" ]
then
	echo "Usage: $0 <geom file> <dir>"
	exit
fi

MESHDIR=$2
VOLGEOM=$1
PARALLEL=YES
T=.commands
if [ "$PARALLEL" == "YES" ]
then
	T=.commands
	rm -f $T
fi

for i in $MESHDIR/*.vtp
do
	B=`basename $i`
	
	CMD="./VTPExtractAll --surf-volgeom=$VOLGEOM $MESHDIR/$B"
	if [ "$PARALLEL" == "YES" ]
	then
		echo $CMD >> $T
	else
		$CMD
	fi
done

if [ "$PARALLEL" == "YES" ]
then
	parallel -j12 --ungroup < $T
	rm -f $T
fi
