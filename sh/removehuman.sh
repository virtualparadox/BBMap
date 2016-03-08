#!/bin/bash
#removehuman in=<infile> out=<outfile>

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified October 17, 2014"
	echo "This script requires at least 11GB RAM."
	echo ""
	echo "Description:  Removes all reads that map to the human genome with at least 95% identity after quality trimming."
	echo ""
	echo "Usage:	removehuman.sh in=<input file> outu=<clean output file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a fasta, fastq, or sam file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "interleaved=auto 	(int) If true, forces fastq input to be paired and interleaved."
	echo "trim=t           	Trim read ends to remove bases with quality below minq."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "untrim=t           	Undo the trimming after mapping."
	echo "minq=4           	Trim quality threshold."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "outm=<file>       	File to output the reads that mapped to human."
	echo ""
	echo "***** All BBMap parameters can be used; run bbmap.sh for more details. *****"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

function removehuman() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA -Xmx10g -cp $CP align2.BBMap minratio=0.9 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/global/projectb/sandbox/gaag/bbtools/hg19 pigz unpigz zl=6 qtrim=rl trimq=10 untrim idtag usemodulo printunmappedcount usejni ztd=2 $@"
	echo $CMD >&2
	$CMD
}

removehuman "$@"
