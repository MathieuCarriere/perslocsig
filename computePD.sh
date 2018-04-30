#!/bin/bash

CATEGORY=airplane
for i in `seq 61 80`;
do
	SHAPE=$CATEGORY/$i.off
	NUM=$(awk NR==2 $SHAPE | sed 's/ .*//')
	echo "Computing persistence diagrams of $SHAPE with $NUM points"

	for j in `seq 0 $((NUM-1))`; do
		PD=OrdPD-$i-$j
		echo "    $PD"
		( ./graph_geodesic $SHAPE $j && ./graph_neighbors $SHAPE ) | ./persUF 0 >> $CATEGORY/$PD

	done
	zip -r $i-OrdPD.zip OrdPD-$i*
	rm OrdPD-$i*
done
