#!/bin/bash

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified March 14, 2014"
	echo ""
	echo "Description:  Generates random synthetic reads from a reference genome.  A read's name indicates its genomic location."
	echo "Allows precise customization of things like insert size and synthetic mutation type, sizes, and rates."
	echo ""
	echo "Usage:	randomreads.sh ref=<reference fasta> out=<output file> minlen=<min length> maxlen=<max length> reads=<number of reads>"
	echo ""
	echo "Optional parameters:"		
	echo "paired=<false>        	Set to true for paired reads."
	echo "interleaved=<false>    	Set to true if paired output is interleaved (rather than in two files)."
	echo "build=<1>             	If multiple references will be used when running in the same working directory, each needs a unique build ID."
	echo "replacenoref=<false>		Set to true to replace N in the reference sequence with random letters."
	echo "seed=<number>			Use this to set the random number generator seed; use -1 for a random seed."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       			This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "					-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
z2="-Xms1g"
EA="-da"
set=0


parseXmx () {
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
			set=1
		elif [[ "$arg" == Xmx* ]]; then
			z="-$arg"
			set=1
		elif [[ "$arg" == -Xms* ]]; then
			z2="$arg"
			set=1
		elif [[ "$arg" == Xms* ]]; then
			z2="-$arg"
			set=1
		elif [[ "$arg" == -da ]] || [[ "$arg" == -ea ]]; then
			EA="$arg"
		fi
	done
}

calcXmx () {
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	
	x=$(ulimit -v)
	#echo "x=$x"
	HOSTNAME=`hostname`
	y=1
	if [[ $x == unlimited ]]; then
		#echo "ram is unlimited"
		echo "This system does not have ulimit set, so max memory cannot be determined.  Attempting to use 4G." 1>&2
		echo "If this fails, please add the argument -Xmx29g (adjusted to ~85 percent of physical RAM)." 1>&2
		y=4
	else
		mult=75;
		if [ $x -ge 1000000000 ]; then
			mult=85
			#echo "ram is 1000g+"
		elif [ $x -ge 500000000 ]; then
			mult=85
			#echo "ram is 500g+"
		elif [ $x -ge 250000000 ]; then
			mult=85
			#echo "ram is 250g+"
		elif [ $x -ge 144000000 ]; then
			mult=85
			#echo "ram is 144g+"
		elif [ $x -ge 120000000 ]; then
			mult=85
			#echo "ram is 120g+"
		elif [ $x -ge 40000000 ]; then
			mult=80
			#echo "ram is 40g+"
		else
			mult=85
			#echo "ram is under 40g"
		fi
		y=$(( ((x-500000)*mult/100)/1000000 ))
	fi
	#echo "y=$y"
	z="-Xmx${y}g"
}
calcXmx "$@"

randomreads() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP align2.RandomReads3 build=1 $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ]; then
	usage
	exit
fi

randomreads "$@"
