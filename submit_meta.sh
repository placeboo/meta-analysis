#!/bin/sh
#
#
cd ~/meta-analysis
for n1 in 500 1000 5000
do
	for n2 in 200 500 1000 5000
	do
		for seed in $(seq 10)
		do
			qsub -cwd ./call_meta.sh $n1 $n2 $seed
		done
	done
done
