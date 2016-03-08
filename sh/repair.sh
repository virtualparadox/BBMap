#!/bin/bash
#repair in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified May 30, 2014"
	echo ""
	echo "Description:  Re-pairs reads that became disordered or had some mates eliminated."
	echo ""
	echo "Usage:	repair.sh in=<input file> out=<pair output file> outs=<singleton output file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "in=<file>        	The 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in."
	echo "in2=<file>       	Use this if 2nd read of pairs are in a different file."
	echo "out=<file>       	The 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out."
	echo "out2=<file>      	Use this to write 2nd read of pairs to a different file."
	echo "outsingle=<file> 	(outs) Write singleton reads here."
	echo ""
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "showspeed=t      	(ss) Set to 'f' to suppress display of processing speed."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "fixinterleaving=f	(fint) Fixes corrupted interleaved files by examining pair names.  Only use on files with broken interleaving."
	echo "repair=t		(rp) Fixes arbitrarily corrupted paired reads by examining read names.  High memory."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx             	This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

repair() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.SplitPairsAndSingles rp $@"
	echo $CMD >&2
	$CMD
}

repair "$@"
