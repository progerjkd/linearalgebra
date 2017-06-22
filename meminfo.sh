#!/bin/bash

if [ ! $1 ]; then
	echo "usage: $0 outputFile"
	exit 1
fi

PID=`pgrep HH`
LOG="$1"
>$LOG

while [ $PID ]; do

	MEM=`pmap -x $PID | grep total | awk '{print $3}'`
#	echo "${MEM} KB" | tee -a $LOG
	MEM=`expr $MEM / 1024`
	echo "${MEM} MB" | tee -a $LOG

	sleep 1

	PID=`pgrep HH`
done

