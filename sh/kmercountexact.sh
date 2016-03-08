#!/bin/bash
#kmercountexact in=<infile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified April 10, 2014"
	echo ""
	echo "Description:  Counts the number of unique kmers in a file."
	echo ""
	echo "Usage:	kmercountexact.sh in=<input>"
	echo ""
	echo "Input may be a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file.  'out' and 'hist' are both optional."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in2=null		Second input file for paired reads"
	echo ""
	echo "Output parameters:"
	echo "khist=<file>        	Print depth histogram"
	echo "out=<file>        	Print kmers and their counts"
	echo "reads=-1		Only process this number of reads, then quit (-1 means all)"
	echo ""
	echo "Hashing parameters:"
	echo "k=31			Kmer length (1 to 31)"
	echo "prefilter=f		If true, use a countmin sketch to store rare kmers more efficiently."
	echo "prehashes=2		Number of hashes for prefilter."
	#TODO echo "minq=6			Ignore kmers containing bases with quality below this."
	#TODO echo "minprob=0.5		Ignore kmers with overall probability of correctness below this."
	echo "threads=X		Spawn exactly X hashing threads (default is number of logical processors).  Total active threads may exceed X due to I/O threads."
	echo ""
	echo "Trimming parameters:"
	echo "ktrim=f          	Trim reads to remove bases matching reference kmers."
	echo "                 	Values: f (don't trim), r (trim right end), l (trim left end), n (convert to N instead of trimming)."
	echo "                 	Any non-whitespace character other than t, f, r, l, n: convert to that symbol rather than trimming, and process short kmers on both ends."
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=4           	Trim quality threshold."
	echo "minavgquality=0  	(maq) Reads with average quality (before trimming) below this will be discarded."
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

kmercountexact() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.KmerCountExact $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

khist "$@"
