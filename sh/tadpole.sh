#!/bin/bash
#tadpole in=<infile>

usage(){
echo "
Written by Brian Bushnell
Last modified May 22, 2015

Description:  Assembles contigs, or extends reads, using short kmers.

Usage:   tadpole.sh in=<file> out=<file>

Input may be fasta or fastq, compressed or uncompressed.

Input parameters:
in=<file>           Primary input file for reads to use as kmer data.
in2=<file>          Second input file for paired data.
extend=<file>       Primary input file for sequences to extend.
extend2=<file>      Second input file for paired reads.
reads=-1            Only process this number of reads, then quit (-1 means all).

Output parameters:
out=<file>          Write contigs (in contig mode).
oute=<file>         Write extended reads (in extend mode).
ihist=<file>        Write insert size histogram (in insert mode).
dump=<file>         Write kmers and their counts.
fastadump=t         Write kmers and counts as fasta versus 2-column tsv.
mincounttodump=1    Only dump kmers with at least this depth.
showstats=t         Print assembly statistics after writing contigs.

Hashing parameters:
k=31                Kmer length (1 to 31)
prealloc=t          Pre-allocate memory rather than dynamically growing; faster and more memory-efficient.  A float fraction (0-1) may be specified, default 1.
prefilter=0         If set to a positive integer, use a countmin sketch to ignore kmers with depth of that value or lower.
prehashes=2         Number of hashes for prefilter.
prefiltersize=0.2   Fraction of memory to use for prefilter.
minq=6              Ignore kmers containing bases with quality below this. (TODO)
minprob=0.5         Ignore kmers with overall probability of correctness below this. (TODO)
threads=X           Spawn X hashing threads (default is number of logical processors).
onepass=f           If true, prefilter will be generated in same pass as kmer counts.  Much faster but counts will be lower, by up to prefilter's depth limit.
rcomp=t             Store and count each kmer together and its reverse-complement.

Assembly parameters:
mincount=3          Minimum kmer depth to assemble.
branchmult1=60      Min ratio of 1st to 2nd-greatest path depth at high depth.
branchmult2=8       Min ratio of 1st to 2nd-greatest path depth at low depth.
branchlower=3       Max value of 2nd-greatest path depth to be considered low.
minextension=2      Do not keep contigs that did not extend at least this much.
mincontig=1         Do not write contigs shorter than this.

Processing modes:
mode=contig         contig: Make contigs from kmers.
                    insert: Measure insert sizes.
                    extend: Extend reads to be longer.

Extension parameters:
extendleft=100      (el) Extend to the left by at most this many bases.
extendright=100     (er) Extend to the right by at most this many bases.
ibb=t               (ignorebackbranches) Do not stop at backward branches.

Quality parameters:
ecc=f               For overlapping paired reads only.  Performs error-
                    correction with BBMerge prior to kmer operations.   

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx15g"
z2="-Xms15g"
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

tadpole() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z $z2 -cp $CP assemble.Tadpole $@"
	echo $CMD >&2
	eval $CMD
}

tadpole "$@"
