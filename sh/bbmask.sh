#!/bin/bash -l
#bbmask in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
calcXmx () {
	x=$(ulimit -v)
	#echo "x=$x"
	HOSTNAME=`hostname`
	y=1
	if [[ $x == unlimited ]] || [[ $HOSTNAME == gpint* ]]; then
		#echo "ram is unlimited"
		echo "This system does not have ulimit set, so max memory cannot be determined.  Attempting to use 1G." 1>&2
		echo "If this fails, please set ulimit or run this program qsubbed or from a qlogin session on Genepool." 1>&2
		y=1000
	fi
	
	mult=30;

	y=$(( ((x-20000)*mult/100)/1000 ))

	if [ $y -ge 2500 ]; then
		y=2500 
	elif [ 200 -ge $y ]; then
		y=200 
	fi
	
	#echo "y=$y"
	z="-Xmx${y}m"
	
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
		fi
	done
}
calcXmx "$@"

bbmask() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java -ea $z -cp $CP jgi.BBDukF k=23 mink=12 ktrim=r qtrim=rl trimq=10 mm=f $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "This script is designed for Genepool nodes."
	echo "Last modified February 25, 2014"
	echo ""
	echo "Description:  Masks sequences of low-complexity, or containing repeat kmers, or covered by mapped reads."
	echo "By default this program will mask using entropy with a window=80 and entropy=0.75"
	echo ""
	echo "Usage:	bbmask.sh in=<file> out=<file> sam=<file,file,...file>"
	echo ""
	echo "Input may be stdin or a fasta or fastq file, raw or gzipped."
	echo "Output may be stdout or a file."
	echo "sam is optional, but may be a comma-delimited list of sam files to mask."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Input sequences to mask. 'in=stdin.fa' will pipe from standard in."
	echo "sam=<file,file>  	Comma-delimited list of sam files.  Optional.  Their mapped coordinates will be masked."
	echo "touppercase=f    	(tuc) Change all letters to upper-case."
	echo "interleaved=auto 	(int) If true, forces fastq input to be paired and interleaved."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo ""
	echo "Output parameters:"
	echo "out=<file>       	Write masked sequences here.  'out=stdout.fa' will pipe to standard out."
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo "lowercase=f        	Set true to convert masked bases to lower case, false to convert them to 'N'."
	echo ""
	echo "Processing parameters:"
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "maskrepeats=f    	Mask areas covered by exact repeat kmers."
	echo "k=5              	Kmer size to use for repeat detection (1-15).  Use mink and maxk to sweep a range of kmers."
	echo "minlen=30     	  	Minimum length of repeat area to mask."
	echo "mincount=3       	Minimum number of repeats to mask."
	echo "masklowentropy=t      Mask areas with low complexity by calculating entropy over a window for a fixed kmer size."
	echo "k2=5     	  	Kmer size to use for entropy calculation (1-15).  Use mink2 and maxk2 to sweep a range.  Large k2 uses more memory."
	echo "window=80     	  	Window size for entropy calculation."
	echo "entropy=0.75      	Mask windows with entropy under this value (0-1).  0 will mask only homopolymers and 1 will mask everything."
	echo "lowercase=f      	Convert masked bases to lower case.  Default is to convert them to N."
	echo ""
	echo ""
	echo "Other parameters:"
	echo "pigz=t    	  	Use pigz to compress.  If argument is a number, that will set the number of pigz threads."
	echo "unpigz=t    	  	Use pigz to decompress."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbmask "$@"
