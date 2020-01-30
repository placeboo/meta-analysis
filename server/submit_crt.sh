#!/bin/sh
#
#
cd ~/meta-analysis
for n in 500 1000 5000
do
	for seed in $(seq 10)
	do
		qsub -cwd ./call_crt.sh $n $seed
	done
done

