#!/bin/bash
#readlength in=<infile>

usage(){
	echo "Generates a length histogram of input reads."
	echo "Written by Brian Bushnell"
	echo "Last modified July 9, 2014"
	echo ""
	echo "Usage:	readlength.sh in=<input file>"
	echo ""
	echo "in=<file>    	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in."
	echo "in2=<file>   	Use this if 2nd read of pairs are in a different file."
	echo "out=<file>   	Write the histogram to this file.  Default is stdout."
	echo "bin=10       	Set the histogram bin size."
	echo "max=80000    	Set the max read length to track."
	echo "round=f      	Places reads in the closest bin, rather than the highest bin of at least readlength."
	echo "nzo=f        	(nonzeroonly) Do not print empty bins."
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

stats() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA -Xmx120m -cp $CP jgi.MakeLengthHistogram $@"
#	echo $CMD >&2
	$CMD
}

stats "$@"
