#!/bin/bash
#bbwrap in=<infile> out=<outfile>

usage(){
	echo "BBWrap v31.x"
	echo "Last modified March 27, 2014"
	echo ""
	echo "Description:  Wrapper for BBMap to allow multiple input and output files for the same reference."
	echo ""
	echo "To index:  	bbwrap.sh ref=<reference fasta>"
	echo "To map:    	bbwrap.sh in=<file,file,...> out=<file,file,...>"
	echo "To map without an index:  	bbwrap.sh ref=<reference fasta> in=in=<file,file,...> out=<file,file,...> nodisk"
	echo ""
	echo "This will probably not work with stdin and stdout."
	echo ""
	echo "Other Parameters:"
	echo ""
	echo "in=<file>        	Input sequences to mask. 'in=stdin.fa' will pipe from standard in."
	echo "mapper=bbmap     	Select mapper.  May be BBMap, BBMapPacBio, or BBMapPacBioSkimmer."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx             		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "                 		-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "This list is not complete.  For more information, please consult $DIR""docs/readme.txt"
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
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
			mult=84
			#echo "ram is 1000g+"
		elif [ $x -ge 500000000 ]; then
			mult=84
			#echo "ram is 500g+"
		elif [ $x -ge 250000000 ]; then
			mult=84
			#echo "ram is 250g+"
		elif [ $x -ge 144000000 ]; then
			mult=84
			#echo "ram is 144g+"
		elif [ $x -ge 120000000 ]; then
			mult=84
			#echo "ram is 120g+"
		elif [ $x -ge 40000000 ]; then
			mult=80
			#echo "ram is 40g+"
		else
			mult=84
			#echo "ram is under 40g"
		fi
		y=$(( ((x-500000)*mult/100)/1000000 ))
	fi
	#echo "y=$y"
	z="-Xmx${y}g"
}
calcXmx "$@"


bbwrap() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP align2.BBWrap build=1 overwrite=true fastareadlen=500 $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbwrap "$@"
