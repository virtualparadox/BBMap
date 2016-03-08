#!/bin/bash
#removesmartbell in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified June 17, 2014"
	echo ""
	echo "Description:  Remove Smart Bell adapters from PacBio reads."
	echo ""
	echo "Usage:	removesmartbell in=<input> out=<output> split=t"
	echo ""
	echo "Input may be a fasta or fastq file, compressed or uncompressed (not H5 files)."
	echo ""
	echo "Parameters:"
	echo "in=file     	Specify the input file, or stdin."
	echo "out=file    	Specify the output file, or stdout."
	echo "split=f     	To split reads at adapters, set split=t.  To mask the adapters, set split=f."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx400m"
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

removesmartbell() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	local CMD="java $EA $z -cp $CP pacbio.RemoveAdapters2 $@"
	echo $CMD >&2
	$CMD
}

removesmartbell "$@"
