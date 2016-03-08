#!/bin/bash
#bbduk in=<infile> out=<outfile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified June 17, 2014"
	echo ""
	echo "Description:  Compares reads to the kmers in a reference dataset, optionally allowing an edit distance."
	echo "Splits the reads into two outputs - those that match the reference, and those that don't."
	echo "Can also trim (remove) the matching parts of the reads rather than binning the reads."
	echo ""
	echo "Usage:	bbduk.sh in=<input file> out=<output file> ref=<contaminant files>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
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
	echo "out=<file>       	(outnonmatch) Write reads here that do not contain kmers matching the database.  'out=stdout.fq' will pipe to standard out."
	echo "out2=<file>      	(outnonmatch2) Use this to write 2nd read of pairs to a different file."
	echo "outmatch=<file>    	(outm or outb) Write reads here that contain kmers matching the database."
	echo "outmatch2=<file>   	(outm2 or outb2) Use this to write 2nd read of pairs to a different file."
	echo "stats=<file>     	Write statistics about which contamininants were detected."
	echo "duk=<file>       	Write statistics in duk's format."
	echo "overwrite=t      	(ow) Grant permission to overwrite files."
	echo "showspeed=t      	(ss) 'f' suppresses display of processing speed."
	echo "ziplevel=2       	(zl) Compression level; 1 (min) through 9 (max)."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	Output quality offset: 33 (Sanger), 64, or auto."
	echo ""
	echo "Histogram output parameters:"
	echo "bhist=<file>		Write a base composition histogram to file."
	echo "qhist=<file>		Write a quality histogram to file."
	echo "indelhist=<file>	Write an indel length histogram to file."
	echo "lhist=<file>		Write a read length histogram to file."
	echo "ehist=<file>		Write an errors-per-read histogram to file."
	echo "gchist=<file>		Write a gc content histogram to file."
	echo "qahist=<file>		Write histogram of match, sub, ins, del by quality score.  Requires sam v1.4 input."
	echo "mhist=<file>		Write histogram of match, substitution, deletion, and insertion rates by read location.  Requires sam v1.4 input."
	echo "idhist=<file>		Write a percent identity histogram to file.  Requires sam v1.4 input."
	echo ""
	echo "Processing parameters:"
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "k=28             	Kmer length used for finding contaminants.  Contaminants shorter than k will not be found.  k must be at least 1."
	echo "rcomp=t            	Look for reverse-complements of kmers in addition to forward kmers."
	echo "maskmiddle=t     	(mm) Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors."
	echo "maxbadkmers=0    	(mbk) Reads with more than this many contaminant kmers will be discarded."
	echo "hammingdistance=0	(hdist) Maximum Hamming distance from ref kmers (subs only).  Memory use is proportional to (3*K)^hdist."
	echo "editdistance=0   	(edist) Maximum edit distance from ref kmers (subs and indels).  Memory use is proportional to (8*K)^edist."
	echo "forbidn=f     	  	(fn) Forbids matching of read kmers containing N.  By default, these will match a reference 'A' if hdist>0 or edist>0, to increase sensitivity."
	echo "minskip=1    	  	(mns) Force minimal skip interval when indexing reference kmers.  1 means use all, 2 means use every other kmer, etc."
	echo "maxskip=99    	  	(mxs) Restrict maximal skip interval when indexing reference kmers."
	echo "                 	Normally all are used for scaffolds<100kb, but with longer scaffolds, up to K-1 are skipped."
	echo "removeifeitherbad=t	(rieb) Paired reads get sent to 'outmatch' if either is match (or either is trimmed shorter than minlen).  Set to false to require both."
	echo "findbestmatch=f  	(fbm) If multiple matches, associate read with sequence sharing most kmers."
	echo ""
	echo "Trimming parameters:"
	echo "ktrim=f          	Trim reads to remove bases matching reference kmers."
	echo "                 	Values: f (don't trim), r (trim right end), l (trim left end), n (convert to N instead of trimming)."
	echo "                 	Any non-whitespace character other than t, f, r, l, n: convert to that symbol rather than trimming, and process short kmers on both ends."
	echo "useshortkmers=f  	(usk) Look for shorter kmers at read tips (only for k-trimming).  Enabling this will disable maskmiddle."
	echo "mink=6           	Minimum length of short kmers.  Setting this automatically sets useshortkmers=t."
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=6           	Trim quality threshold."
	echo "minlength=20     	(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter."
	echo "minavgquality=0  	(maq) Reads with average quality (before trimming) below this will be discarded."
	echo "otm=f            	(outputtrimmedtomatch) Output reads trimmed to shorter than minlength to outm rather than discarding."
	echo "tp=0             	(trimpad) Trim this much extra around matching kmers."
	echo "tbo=f            	(trimbyoverlap) Trim adapters based on where paired reads overlap."
	echo "tpe=f            	(trimpairsevenly) When kmer right-trimming, trim both reads to the minimum length of either."
	echo "forcetrimleft=0     	(ftl) If nonzero, trim left bases of the read to this position (exclusive, 0-based)."
	echo "forcetrimright=0    	(ftr) If nonzero, trim right bases of the read after this position (exclusive, 0-based)."
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
	freeRam 2000m 42
	z="-Xmx${RAM}m"
}
calcXmx "$@"

bbduk() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	#module load samtools
	local CMD="java $EA $z -cp $CP jgi.BBDukF $@"
	echo $CMD >&2
	$CMD
}

bbduk "$@"
