#!/bin/bash
#mapPacBio in=<infile> out=<outfile> ref=<reference>

usage(){
	bash "$DIR"bbmap.sh
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

mapPacBio() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java -ea $z -cp $CP align2.BBMapPacBio build=1 overwrite=true minratio=0.40 fastareadlen=500 ambiguous=best minscaf=100 startpad=4000 stoppad=4000 midpad=1000 $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ]; then
	usage
	exit
fi

mapPacBio "$@"
