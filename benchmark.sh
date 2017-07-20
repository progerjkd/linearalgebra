#!/bin/bash

# run as: ./benchmark.sh | tee log
rm input/*\.memory

for INPUT in `ls -Sr input/*_sym* | grep -v memory`; do
	echo -e "\nInput file: ${INPUT}"
	time ./runHHSimetrica2 ${INPUT} &
	./meminfo.sh "${INPUT}.memory" >/dev/null &
	PID=`pgrep HH`
	while [ $PID ]; do
		sleep 1
		PID=`pgrep HH`
	done
	sleep 10


done

for INPUT in `ls -Sr input/*_asym* | grep -v memory`; do
	echo -e "\nInput file: ${INPUT}"
	time ./runHHAssimetrica2 ${INPUT}
	./meminfo.sh "${INPUT}.memory" >/dev/null &
	PID=`pgrep HH`
	while [ $PID ]; do
		sleep 1
		PID=`pgrep HH`
	done
	sleep 10
done

