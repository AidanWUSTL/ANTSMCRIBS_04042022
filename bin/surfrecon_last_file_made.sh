#!/bin/bash

for i in SurfReconDeformable/*
do
	B=`basename $i`
	W=`ls -1tr $i/temp/*.vtp | tail -n 1`
	D=`basename $W`
	echo -e "$B\t$D"
done
