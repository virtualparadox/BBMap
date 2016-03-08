#!/bin/bash
#translate6frames in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified April 3, 2014"
	echo ""
	echo "Description:  Translates nucleotide sequences to all 6 amino acid frames."
	echo ""
	echo "Usage:	translate6frames.sh in=<input file> out=<output file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Main input. in=stdin.fq will pipe from stdin."
	echo "in2=<file>       	Input for 2nd read of pairs in a different file."
	echo "interleaved=auto 	(int) t/f overrides interleaved autodetection."
	echo "qin=auto         	Input quality offset: 33 (Sanger), 64, or auto."
	echo "reads=-1         	If positive, quit after processing X reads or pairs."
	echo ""
	echo "Output parameters:"
	echo "out=<file>       	(outnonmatch) Write reads here that do not contain kmers matching the database.  'out=stdout.fq' will pipe to standard out."
	echo "out2=<file>      	(outnonmatch2) Use this to write 2nd read of pairs to a different file."
	echo "overwrite=t      	(ow) Grant permission to overwrite files."
	echo "ziplevel=2       	(zl) Compression level; 1 (min) through 9 (max)."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	Output quality offset: 33 (Sanger), 64, or auto."
	echo "tag=t            	Tag read id with the frame, adding e.g. ' fr1'"
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
		echo "This system does not have ulimit set, so max memory cannot be determined.  Attempting to use 2G." 1>&2
		echo "If this fails, please add the argument -Xmx29g (adjusted to ~85 percent of physical RAM)." 1>&2
		y=2000
	fi
	
	mult=42;

	y=$(( ((x-20000)*mult/100)/1000 ))

	if [ $y -ge 15000 ]; then
		y=15000
	elif [ 200 -ge $y ]; then
		y=200
	fi
	
	#echo "y=$y"
	z="-Xmx${y}m"
}
calcXmx "$@"

translate6frames() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP jgi.TranslateSixFrames $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ]; then
	usage
	exit
fi

translate6frames "$@"
