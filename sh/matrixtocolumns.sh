#!/bin/bash
#matrixtocolumns in1=<infile> in2=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified November 23, 2014"
	echo ""
	echo "Description:  Turns identity matrices into 2-column format for plotting."
	echo ""
	echo "Usage:	matrixtocolumns.sh in1=<matrix1> in2=<matrix2> out=<file>"
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx2g"
z2="-Xms2g"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

matrixtocolumns() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP driver.CorrelateIdentity $@"
	echo $CMD >&2
	$CMD
}

matrixtocolumns "$@"
