#!/bin/bash
#merge in=<infile> out=<outfile>

function usage(){
	echo "BBMerge v4.2"
	echo "Written by Brian Bushnell"
	echo "Last modified June 24, 2014"
	echo ""
	echo "Description:  Merges paired reads into single reads by overlap detection."
	echo "With sufficient coverage, can also merge nonoverlapping reads using gapped kmers."
	echo ""
	echo "Usage for interleaved files:	bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>"
	echo "Usage for paired files:     	bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or scarf file, raw or gzipped."
	echo "Output may be stdout or a file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=null  		Primary input. 'in2' will specify a second file."
	echo "interleaved=auto  	May be set to true or false to override autodetection of whether the input file as interleaved."
	echo "reads=-1  		Only process this number of read pairs, then quit (-1 means all)."
	echo ""
	echo "Output parameters:"
	echo "out=<file>        	File for merged reads. 'out2' will specify a second file."
	echo "outu=<file>       	File for unmerged reads. 'outu2' will specify a second file."
	echo "outinsert=<file>    	File list of read names and their insert sizes."
	echo "hist=null    		Insert length histogram output file."
	echo "ziplevel=2   		Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo ""
	echo "Trimming parameters:"
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq.  Performed BEFORE merging."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=10           	Trim quality threshold."
	echo "minlength=20     	(ml) Reads shorter than this after trimming (before merging) will be discarded."
	echo "           		Pairs will be discarded only if both are shorter."
	echo "trimonfailure=t     	(tof) If detecting insert size by overlap fails, the reads will be trimmed and this will be re-attempted."
	echo "tbo=f            	(trimbyoverlap) Trim overlapping reads to remove rightmost (3') non-overlapping portion, instead of joining."
	echo ""
	echo "Processing Parameters:"
	echo "join=t			Create merged reads.  If set to false, you can still generate an insert histogram."
	echo "useoverlap=t		Attempt merge based on paired read overlap."
	echo "minoverlap=12		Minimum number of overlapping bases to consider merging reads."
	echo "minoverlap0=8		Minimum number of overlapping bases to consider for deciding ambiguity."
	echo "minoverlapinsert=25	Do not look for insert sizes shorter than this."
	echo "mininsert=35		Reads with insert sizes less than this (after merging) will be discarded."
	echo "minqo=9    		Ignore bases with quality below this."
	echo "margin=2    		The best overlap must have at least 'margin' fewer mismatches than the second best."
	echo ""
	echo "Processing Modes (these are mutually exclusive macros that set other parameters):"
	echo "strict=f    		Decrease false positive rate and merging rate."
	echo "fast=f      		Increase speed and slightly decrease merging rate."
	echo "normal=t    		Default."
	echo "loose=f     		Increase false positive rate and merging rate."
	echo "vloose=f    		Greatly increase false positive rate and merging rate."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "           		For example, -Xmx400m will specify 400 MB RAM."
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

function merge() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.BBMerge $@"
	echo $CMD >&2
	$CMD
}

merge "$@"
