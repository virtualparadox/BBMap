#!/bin/bash
#testformat in=<infile>

usage(){
	echo "Tests the format of a sequence file based on name and contents."
	echo "Written by Brian Bushnell"
	echo "Last modified October 23, 2014"
	echo ""
	echo "Description:  Tests the file extensions and contente of filee to determine format, quality, compression, and interleaving."
	echo ""
	echo "Usage:	testformat.sh <file>"
	echo ""
	echo "More than one file may be specified."
	echo "Note that ASCII-33 (sanger) and ASCII-64 (illumina) cannot always be differentiated."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx40m"
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

testformat() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA $z -cp $CP fileIO.FileFormat $@"
#	echo $CMD >&2
	$CMD
}

testformat "$@"
