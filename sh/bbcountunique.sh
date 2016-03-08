#!/bin/bash
#bbuniqueness in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified March 24, 2014"
	echo ""
	echo "Description:  Generates a kmer uniqueness histogram, binned by file position."
	echo "There are 3 columns for single reads, 6 columns for paired:"
	echo "count      	number of reads or pairs processed"
	echo "r1_first   	percent unique 1st kmer of read 1"
	echo "r1_rand    	percent unique random kmer of read 1"
	echo "r2_first   	percent unique 1st kmer of read 2"
	echo "r2_rand    	percent unique random kmer of read 2"
	echo "pair       	percent unique concatenated kmer from read 1 and 2"
	echo ""
	echo ""
	echo "Usage:	bbuniqueness.sh in=<input> out=<output>"
	echo ""
	echo "Input may be a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in2=null		Second input file for paired reads"
	echo "interleaved=auto	May be set to true or false to force the input read file to ovverride autodetection of the input file as paired interleaved."
	echo "samplerate=1		Set to below 1 to sample a fraction of input reads."
	echo "reads=-1		Only process this number of reads, then quit (-1 means all)"
	echo ""
	echo "Output parameters:"
	echo "out=<file>        	File for output stats"
	echo ""
	echo "Hashing parameters:"
	echo "k=20			Kmer length (values under 32 are most efficient, but arbitrarily high values are supported)"
	echo ""
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

z="-Xmx1g"
z2="-Xms1g"
EA="-ea"
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

bbuniqueness() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.CalcUniqueness $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbuniqueness "$@"
