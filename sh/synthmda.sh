#!/bin/bash
#synthmda in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified November 2, 2014"
	echo ""
	echo "Description:  Generates synthetic reads following an MDA-amplified singe cell's coverage distribution."
	echo ""
	echo "Usage:	synthmda.sh in=<reference> out=<reads out file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "reads=12000000 	Generate this many reads."
	echo "paired=t    		Generate paired reads."
	echo "length=150  		Reads should be this long."
	echo "minlen=10000		Min length of MDA contig."
	echo "maxlen=150000		Max length of MDA contig."
	echo "cycles=9    		Number of MDA cycles (higher is more spiky)."
	echo "initialratio=1.3	Fraction of genome initially replicated (lower is more spiky)."
	echo "ratio=1.7   		Fraction of genome replicated per cycle."
	echo "refout=null 		Write MDA'd genome to this file."
	echo "perfect=0   		This fraction of reads will be error-free."
	echo "amp=200     		'amp' flag sent to RandomReads (higher is more spiky)."
	echo "build=7     		Index MDA'd genome in this location."
	echo "prefix=null 		Generated reads will start with this prefix."
	echo "overwrite=t 		(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "           		-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 4000m 80
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

synthmda() {
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.SynthMDA $@"
	echo $CMD >&2
	$CMD
}

synthmda "$@"
