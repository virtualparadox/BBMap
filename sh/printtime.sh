#!/bin/bash
#printtime in=<infile> out=<outfile>

function usage(){
	echo "Prints time elapsed since last called on the same file."
	echo "Written by Brian Bushnell"
	echo "Last modified May 23, 2014"
	echo ""
	echo "Usage:	printtime.sh <filename>"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

function printtime() {
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA -Xmx8m -cp $CP align2.PrintTime $@"
	echo $CMD >&2
	$CMD
}

printtime "$@"
