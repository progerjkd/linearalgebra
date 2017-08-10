#!/bin/bash

# run as: ./benchmark.sh | tee log
rm input/*\.memory

for INPUT in `ls -Sr input/*_sym* | grep -v memory`; do
	echo -e "\nInput file: ${INPUT}"
	time ./runHHSimetrica3 ${INPUT} &
	./meminfo.sh "${INPUT}.memory" >/dev/null &
	PID=`pgrep HH`
	while [ $PID ]; do
		sleep 1
		PID=`pgrep HH`
	done
	echo -e "\nInput file: ${INPUT}"
	sleep 1


done

for INPUT in `ls -Sr input/*_asym* | grep -v memory`; do
	echo -e "\nInput file: ${INPUT}"
	time ./runHHAssimetrica3 ${INPUT}
	./meminfo.sh "${INPUT}.memory" >/dev/null &
	PID=`pgrep HH`
	while [ $PID ]; do
		sleep 1
		PID=`pgrep HH`
	done
	echo -e "\nInput file: ${INPUT}"
	sleep 1
done

