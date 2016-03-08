#!/bin/bash
#merge in=<infile> out=<outfile>

function usage(){
echo "
BBMerge v7.3
Written by Brian Bushnell and Jonathan Rood
Last modified May 27, 2015

Description:  Merges paired reads into single reads by overlap detection.
With sufficient coverage, can also merge nonoverlapping reads using gapped kmers.

Usage for interleaved files:	bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
Usage for paired files:     	bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

Input may be stdin or a fasta, fastq, or scarf file, raw or gzipped.
Output may be stdout or a file.


Optional parameters (and their defaults)

Input parameters:
in=null              Primary input. 'in2' will specify a second file.
interleaved=auto     May be set to true or false to override autodetection of
                     whether the input file as interleaved.
reads=-1             Quit after this many read pairs (-1 means all).


Output parameters:
out=<file>           File for merged reads. 'out2' will specify a second file.
outu=<file>          File for unmerged reads. 'outu2' will specify a second file.
outadapter=<file>    (outa) File to write consensus adapter sequences.
outinsert=<file>     File list of read names and their insert sizes.
hist=null            Insert length histogram output file.
nzo=t                Only print histogram bins with nonzero values.
ziplevel=2           Set to 1 (lowest) through 9 (max) to change compression
                     level; lower compression is faster.
ordered=f            Output reads in same order as input.
mix=f                Output both the merged (or mergable) and unmerged reads
                     in the same file (out=).  Useful for ecc mode.


Trimming/Filtering parameters:
qtrim=f              Trim read ends to remove bases with quality below minq.  
                     Performed BEFORE merging.
                     Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=10             Trim quality threshold.
minlength=0          (ml) Reads shorter than this after trimming, but before
                     merging, will be discarded. Pairs will be discarded only
                     if both are shorter.
maxlength=-1         Reads with longer insert sizes will be discarded.
trimonfailure=t      (tof) If detecting insert size by overlap fails,
                     the reads will be trimmed and this will be re-attempted.
tbo=f                (trimbyoverlap) Trim overlapping reads to remove 
                     rightmost (3') non-overlapping portion, instead of joining.
minavgquality=0      (maq) Reads with average quality below this (AFTER trimming)
                     will not be attempted to be merged.
maxexpectederrors=0  (mee) If positive, reads with more combined expected 
                     errors than this will not be attempted to be merged.


Processing Parameters:
usejni=f             (jni) Do overlapping in C code, which is faster.  Requires
                     compiling the C code; details are in /jni/README.txt.
merge=t              Create merged reads.  If set to false, you can still 
                     generate an insert histogram.
ecc=f                Error-correct the overlapping part, but don't merge.
useoverlap=t         Attempt find the insert size using read overlap.
mininsert=35         Minimum insert size to merge reads.
mininsert0=35        Insert sizes less than this will not be considered.
                     Must be less than or equal to mininsert.
minoverlap=12        Minimum number of overlapping bases to allow merging.
minoverlap0=8        Overlaps shorter than this will not be considered.
                     Must be less than or equal to minoverlap.
minq=9               Ignore bases with quality below this.
maxq=41              Cap output quality scores at this.
entropy=t            Increase the minimum overlap requirement for low-
                     complexity reads.
efilter=6            Ban overlaps with over this many times the expected 
                     number of errors.  Lower is more strict.
pfilter=0.00002      Ban improbable overlaps.  Higher is more strict. 0 will
                     disable the filter; 1 will allow only perfect overlaps.
kfilter=f            Ban overlaps that create novel kmers.  Does not help.
lowercase=f          Expect lowercase letters to signify adapter sequence.
ouq=f                Calculate best overlap using quality values.
owq=t                Calculate best overlap without using quality values.
usequality=t         If disabled, quality values are completely ignored,
                     both for overlap detection and filtering.  May be useful
                     for data with inaccurate quality values.
iupacton=f           (itn) Change ambiguous IUPAC symbols to N.


Normal Mode: 
normalmode=f         Original BBMerge algorithm.  Faster, but lower overall
                     merge rate.
margin=2             The best overlap must have at least 'margin' fewer 
                     mismatches than the second best.
mismatches=3         Allow up to this many mismatches in the overlap.
requireratiomatch=f  (rrm) Require the answer from normal mode and ratio mode
                     to agree, reducing false positives if both are enabled.


Ratio Mode: 
ratiomode=t          Newer algorithm.  Slower, but higher merge rate.
                     Much better for long overlaps and high error rates.
maxratio=0.09        Max error rate; higher increases merge rate.
ratiomargin=5.5      Lower increases merge rate; min is 1.
ratiooffset=0.55     Lower increases merge rate; min is 0.
ratiominoverlapreduction=3  This is the difference between minoverlap in 
                     normal mode and minoverlap in ratio mode; generally, 
                     minoverlap should be lower in ratio mode.

*** Ratio Mode and Normal Mode may be used alone or simultaneously ***


Strictness (these are mutually exclusive macros that set other parameters):
strict=f             Decrease false positive rate and merging rate.
verystrict=f         (vstrict) Greatly decrease FP and merging rate.
ultrastrict=f        (ustrict) Decrease FP and merging rate even more.
maxstrict=f          (xstrict) Maximally decrease FP and merging rate.
loose=f              Increase false positive rate and merging rate.
veryloose=f          (vloose) Greatly increase FP and merging rate.
ultraloose=f         (uloose) Increase FP and merging rate even more.
maxloose=f           (xloose) Maximally decrease FP and merging rate.
fast=f               Fastest possible mode; less accurate.


Java Parameters:
-Xmx                 This will be passed to Java to set memory usage, 
                     overriding the program's automatic memory detection.
                     For example, -Xmx400m will specify 400 MB RAM.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

z="-Xmx1000m"
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
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP jgi.BBMerge $@"
	echo $CMD >&2
	eval $CMD
}

merge "$@"
