#!/bin/bash

# run as follow:
# ./benchmark.sh 2>&1 |tee log

for INPUT in `ls -Sr input/*_sym*`; do
	echo -e "\nInput file: ${INPUT}"
	time ./runHHSimetrica2 ${INPUT}
done

for INPUT in `ls -Sr input/*_asym*`; do
	echo -e "\nInput file: ${INPUT}"
	time ./runHHAssimetrica2 ${INPUT}
done

