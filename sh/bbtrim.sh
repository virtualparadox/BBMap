#!/bin/bash
#bbtrim in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified March 14, 2014"
	echo ""
	echo "Description:  Performs quality-trimming and/or kmer-trimming on reads."
	echo ""
	echo "Usage:	bbtrim.sh in=<input file> out=<output file> ref=<adapter files> trimq=10 qtrim=rl"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo "ref is optional, but may be a comma-delimited list of adapter/primer/linker fasta files."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	Main input. in=stdin.fq will pipe from stdin."
	echo "in2=<file>       	Input for 2nd read of pairs in a different file."
	echo "ref=<file,file>  	Comma-delimited list of reference files."
	echo "touppercase=f    	(tuc) Change all bases upper-case."
	echo "interleaved=auto 	(int) t/f overrides interleaved autodetection."
	echo "qin=auto         	Input quality offset: 33 (Sanger), 64, or auto."
	echo "reads=-1         	If positive, quit after processing X reads or pairs."
#	echo "skipreads=0      	Ignore this many initial reads (or pairs) and process the rest."
	echo ""
	echo "Output parameters:"
	echo "out=<file>       	Output file.  out=stdout.fq will pipe to stdout."
	echo "out2=<file>      	Output for 2nd read of pairs."
	echo "outbad=<file>    	(outb) Output for reads shorter than minlen."
	echo "outbad2=<file>   	(outb2) Output for 2nd read of pairs."
	echo "overwrite=t      	(ow) Grant permission to overwrite files."
	echo "showspeed=t      	(ss) 'f' suppresses display of processing speed."
	echo "ziplevel=2       	(zl) Compression level; 1 (min) through 9 (max)."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	Output quality offset: 33 (Sanger), 64, or auto."
	echo ""
	echo "Processing parameters:"
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "rcomp=t            	Look for reverse-complements of kmers in addition to forward kmers."
	echo "hammingdistance=0	(hdist) Maximum Hamming distance from ref kmers (subs only).  Memory use is proportional to (3*K)^hdist."
	echo "editdistance=0   	(edist) Maximum edit distance from ref kmers (subs and indels).  Memory use is proportional to (8*K)^edist."
	echo "forbidn=f     	  	(fn) Forbids matching of read kmers containing N.  By default, these will match a reference 'A' if hdist>0 or edist>0, to increase sensitivity."
	echo "minskip=1    	  	(mns) Force minimal skip interval when indexing reference kmers.  1 means use all, 2 means use every other kmer, etc."
	echo "maxskip=99    	  	(mxs) Restrict maximal skip interval when indexing reference kmers."
	echo "                 	Normally all are used for scaffolds<100kb, but with longer scaffolds, up to K-1 are skipped."
	echo "removeifeitherbad=t	(rieb) Paired reads get sent to 'outmatch' if either is match (or either is trimmed shorter than minlen).  Set to false to require both."
	echo ""
	echo "Trimming parameters:"
	echo "k=25             	Kmer length used for finding contaminants.  Contaminants shorter than k will not be found.  k must be at least 1."
	echo "ktrim=r          	Trim reads to remove bases matching reference kmers."
	echo "                 	Values: f (filter, don't trim), r (trim right end), l (trim left end), n (convert to N instead of trimming)."
	echo "                 	Any single non-whitespace character other than t, f, r, l, n: convert to that symbol rather than trimming."
	echo "useshortkmers=t  	(usk) Look for shorter kmers at read tips (only for k-trimming).  Enabling this will disable maskmiddle."
	echo "mink=12           	Minimum length of short kmers.  Setting this automatically sets useshortkmers=t."
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=10           	Trim quality threshold."
	echo "minlength=20     	(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter."
	echo "minavgquality=0  	(maq) Reads with average quality (before trimming) below this will be discarded."
	echo "otm=f            	(outputtrimmedtomatch) Output reads trimmed to shorter than minlength to outm rather than discarding."
	echo ""
#	echo "Other parameters:"
#	echo "array=t    	  	Use HashArray data structure."
#	echo "forest=f    	  	Use HashForest data structure."
#	echo "table=f    	  	Use KmerTable data structure."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "There is a changelog at /global/projectb/sandbox/gaag/bbtools/docs/changelog_bbduk.txt"
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
z2="-Xms200m"
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
		echo "This system does not have ulimit set, so max memory cannot be determined.  Attempting to use 1G." 1>&2
		echo "If this fails, please add the argument -Xmx29g (adjusted to ~85 percent of physical RAM)." 1>&2
		y=1000
	fi
	
	mult=30;

	y=$(( ((x-20000)*mult/100)/1000 ))

	if [ $y -ge 1500 ]; then
		y=1500 
	elif [ 200 -ge $y ]; then
		y=200 
	fi
	
	#echo "y=$y"
	z="-Xmx${y}m"
}
calcXmx "$@"

bbduk() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP jgi.BBDukF k=23 mink=12 ktrim=r qtrim=rl trimq=10 mm=f $@"
	echo $CMD >&2
	$CMD
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbduk "$@"
