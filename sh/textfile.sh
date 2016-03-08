#!/bin/bash

function usage(){
	echo "Displays contents of a file."
	echo "Written by Brian Bushnell"
	echo "Last modified May 23, 2014"
	echo ""
	echo "Usage:	textfile.sh <file> <start line> <stop line>"
	echo ""
	echo "Start line and stop line are zero-based, and may be omitted."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

function tf() {
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA -Xmx120m -cp $CP fileIO.TextFile $@"
	echo $CMD >&2
	$CMD
}

tf "$@"
