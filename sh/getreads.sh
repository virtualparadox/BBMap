#!/bin/bash

function usage(){
	echo "Selects reads with designated numeric IDs."
	echo "Written by Brian Bushnell"
	echo "Last modified June 26, 2014"
	echo ""
	echo "Usage:	getreads.sh in=<file> id=<number,number,number...> out=<file>"
	echo ""
	echo "The first read (or pair) has ID 0, the second read (or pair) has ID 1, etc."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

function tf() {
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA $z -cp $CP jgi.GetReads $@"
	echo $CMD >&2
	$CMD
}

tf "$@"