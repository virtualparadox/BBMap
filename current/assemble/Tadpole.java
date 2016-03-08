package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicLong;

import jgi.BBMerge;

import bloom.KCountArray;
import bloom.KmerCount7MTA;
import bloom.KmerCountAbstract;

import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.HashBuffer;
import kmer.Primes;

import stream.ByteBuilder;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;

import align2.ListNum;
import align2.LongList;
import align2.ReadComparatorID;
import align2.ReadLengthComparator;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.CoverageArray3;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;


/**
 * Short-kmer assembler based on KmerCountExact.
 * @author Brian Bushnell
 * @date May 15, 2015
 *
 */
public class Tadpole {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		Tadpole cke=new Tadpole(args);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		cke.process(t);
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		outstream.println("Syntax:\nTODO"); //TODO
//		outstream.println("\njava -ea -Xmx20g -cp <path> jgi.BBAsm in=<input file>");
//		outstream.println("\nOptional flags:");
//		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
//		outstream.println("in2=<file>         \tUse this if 2nd read of pairs are in a different file.");
//		outstream.println("out=<file>         \tDump kmers and counts to this file.");
//		outstream.println("");
//		outstream.println("threads=auto       \t(t) Set number of threads to use; default is number of logical processors.");
//		outstream.println("overwrite=t        \t(ow) Set to false to force the program to abort rather than overwrite an existing file.");
//		outstream.println("showspeed=t        \t(ss) Set to 'f' to suppress display of processing speed.");
//		outstream.println("interleaved=auto   \t(int) If true, forces fastq input to be paired and interleaved.");
//		outstream.println("k=31               \tKmer length used for finding contaminants.  Contaminants shorter than k will not be found.");
//		outstream.println("minavgquality=0    \t(maq) Reads with average quality (before trimming) below this will be discarded.");
//		outstream.println("touppercase=f      \t(tuc) Change all letters in reads and reference to upper-case.");
//		outstream.println("qtrim=f            \tTrim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers. ");
//		outstream.println("                   \tValues: t (trim both ends), f (neither end), r (right end only), l (left end only).");
//		outstream.println("minq=4             \tTrim quality threshold.");
//		outstream.println("minlength=2        \t(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.");
//		outstream.println("ziplevel=2         \t(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.");
//		outstream.println("fastawrap=70       \tLength of lines in fasta output");
//		outstream.println("qin=auto           \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
//		outstream.println("qout=auto          \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
//		outstream.println("rcomp=t            \tLook for reverse-complements of kmers also.");
//		outstream.println("forest=t           \tUse HashForest data structure");
//		outstream.println("table=f            \tUse KmerTable data structure");
//		outstream.println("array=f            \tUse HashArray data structure");
	}
	
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Tadpole(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		ByteFile.FORCE_MODE_BF2=Shared.threads()>2;
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		boolean rcomp_=true;
		boolean useForest_=false, useTable_=false, useArray_=true, prealloc_=false;
		long skipreads_=0;
		int k_=31;
		int ways_=-1;
		int filterMax_=2;
		boolean ecc_=false;
		boolean useOwnership_=true;
		
		{
			boolean b=false;
			assert(b=true);
			EA=b;
		}
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
//			}else if(parser.parseTrim(arg, a, b)){
//				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in1.add(ss);
					}
				}
			}else if(a.equals("in2")){
				in2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in2.add(ss);
					}
				}
			}else if(a.equals("ine") || a.equals("ine1") || a.equals("extend") || a.equals("extend1")){
				ine1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						ine1.add(ss);
					}
				}
			}else if(a.equals("ine2") || a.equals("extend2")){
				ine2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						ine2.add(ss);
					}
				}
			}else if(a.equals("oute") || a.equals("oute1") || a.equals("extend") || a.equals("extend1")){
				oute1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						oute1.add(ss);
					}
				}
			}else if(a.equals("oute2") || a.equals("extend2")){
				oute2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						oute2.add(ss);
					}
				}
			}else if(a.equals("out") || a.equals("contigs")){
				outContigs=b;
			}else if(a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
			}else if(a.equals("mincounttodump")){
				minToDump=(int)Tools.parseKMG(b);
			}else if(a.equals("hist") || a.equals("khist")){
				outHist=b;
			}else if(a.equals("ihist") || a.equals("inserthistogram")){
				outInsert=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=(int)Tools.parseKMG(b);
			}else if(a.equals("mode")){
				if(Character.isDigit(b.charAt(0))){
					processingMode=(int)Tools.parseKMG(b);
				}else if(b.equalsIgnoreCase("contig")){
					processingMode=contigMode;
				}else if(b.equalsIgnoreCase("extend")){
					processingMode=extendMode;
				}else if(b.equalsIgnoreCase("insert")){
					processingMode=insertMode;
				}else{
					assert(false) : "Unknown mode "+b;
				}
			}else if(a.equals("ownership")){
				useOwnership_=Tools.parseBoolean(b);
			}else if(a.equals("showstats") || a.equals("stats")){
				showStats=Tools.parseBoolean(b);
			}else if(a.equals("maxextension") || a.equals("maxe")){
				extendLeft=extendRight=(int)Tools.parseKMG(b);
			}else if(a.equals("extendright") || a.equals("er")){
				extendRight=(int)Tools.parseKMG(b);
			}else if(a.equals("extendleft") || a.equals("el")){
				extendLeft=(int)Tools.parseKMG(b);
			}else if(a.equals("minextension") || a.equals("mine")){
				minExtension=(int)Tools.parseKMG(b);
			}else if(a.equals("maxcontiglength") || a.equals("maxcontig") || a.equals("maxlength") || a.equals("maxlen") || a.equals("maxc")){
				maxContigLen=(int)Tools.parseKMG(b);
				if(maxContigLen<0){maxContigLen=1000000000;}
			}else if(a.equals("mincontiglength") || a.equals("mincontig") || a.equals("minlength") || a.equals("minlen") || a.equals("minc")){
				minContigLen=(int)Tools.parseKMG(b);
			}else if(a.equals("branchlower") || a.equals("branchlowerconst")){
				branchLowerConst=(int)Tools.parseKMG(b);
			}else if(a.equals("branchmult2")){
				branchMult2=(int)Tools.parseKMG(b);
			}else if(a.equals("branchmult1")){
				branchMult1=(int)Tools.parseKMG(b);
			}else if(a.equals("mincount") || a.equals("mincov") || a.equals("mindepth")){
				minCount=(int)Tools.parseKMG(b);
			}else if(a.equals("forest")){
				useForest_=Tools.parseBoolean(b);
				if(useForest_){useTable_=useArray_=false;}
			}else if(a.equals("table")){
				useTable_=Tools.parseBoolean(b);
				if(useTable_){useForest_=useArray_=false;}
			}else if(a.equals("array")){
				useArray_=Tools.parseBoolean(b);
				if(useArray_){useTable_=useForest_=false;}
			}else if(a.equals("ways")){
				ways_=(int)Tools.parseKMG(b);
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
				buflen=(int)Tools.parseKMG(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nk needs an integer value from 1 to 31, such as k=27.  Default is 31.\n";
				k_=(int)Tools.parseKMG(b);
			}else if(a.equals("skipreads")){
				skipreads_=Tools.parseKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("ecc")){
				ecc_=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				verbose=Tools.parseBoolean(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isAlphabetic(b.charAt(0))){
					prealloc_=Tools.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc_=(preallocFraction>0);
				}
			}else if(a.equals("prefilter")){
				if(b==null || b.length()<1 || !Character.isDigit(b.charAt(0))){
					prefilter=Tools.parseBoolean(b);
				}else{
					filterMax_=(int)Tools.parseKMG(b);
					prefilter=filterMax_>0;
				}
			}else if(a.equals("prefiltersize")){
				prefilterFraction=Tools.max(0, Double.parseDouble(b));
				assert(prefilterFraction<=1) : "prefiltersize must be 0-1, a fraction of total memory.";
				prefilter=prefilterFraction>0;
			}else if(a.equals("prehashes") || a.equals("hashes")){
				prehashes=(int)Tools.parseKMG(b);
			}else if(a.equals("onepass")){
				onePass=Tools.parseBoolean(b);
			}else if(a.equals("passes")){
				int passes=(int)Tools.parseKMG(b);
				onePass=(passes<2);
			}else if(a.equals("histcolumns")){
				histColumns=(int)Tools.parseKMG(b);
			}else if(a.equals("histmax")){
				histMax=(int)Tools.parseKMG(b);
			}else if(a.equals("histheader")){
				histHeader=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				histZeros=!Tools.parseBoolean(b);
			}else if(a.equals("ilb") || a.equals("ignoreleftbranches") || a.equals("ignoreleftjunctions") || a.equals("ibb") || a.equals("ignorebackbranches")){
				extendThroughLeftJunctions=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(processingMode==extendMode){
			if(extendLeft==-1){extendLeft=100;}
			if(extendRight==-1){extendRight=100;}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
//			qtrimLeft=parser.qtrimLeft;
//			qtrimRight=parser.qtrimRight;
//			trimq=parser.trimq;
//			
//			minAvgQuality=parser.minAvgQuality;
//			minAvgQualityBases=parser.minAvgQualityBases;
		}
		
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
//			long tmemory=Runtime.getRuntime().totalMemory();
			usableMemory=(long)Tools.max(((memory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.75)), memory*0.45);
			filterMemory=(long)(prefilter ? usableMemory*prefilterFraction : 0);
			tableMemory=(long)(usableMemory*.95-filterMemory);
		}
		
		if(ways_<1){
			long maxKmers=(2*tableMemory)/12;
			long minWays=Tools.min(10000, maxKmers/Integer.MAX_VALUE);
			ways_=(int)Tools.max(31, THREADS*2, minWays);
			ways_=(int)Primes.primeAtLeast(ways_);
			assert(ways_>0);
			System.err.println("ways="+ways_);
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		onePass=onePass&prefilter;
		prealloc=prealloc_;
		useForest=useForest_;
		useTable=useTable_;
		useArray=useArray_;
		rcomp=rcomp_;
		WAYS=ways_;
		filterMax=Tools.min(filterMax_, 0x7FFFFFFF);
		ecc=ecc_;
		useOwnership=useOwnership_;
		
		k=k_;
		k2=k-1;
		
		if(k<1 || k>31){throw new RuntimeException("\nk needs an integer value from 1 to 31, such as k=27.  Default is 31.\n");}
		
		if(initialSize<1){
			int mult=12;
			if(useOwnership){mult+=4;}
			if(processingMode==extendMode){mult+=1;}
			else if(processingMode==contigMode){mult+=1;}
			
			final long memOverWays=tableMemory/(mult*WAYS);
			final double mem2=(prealloc ? preallocFraction : 1)*tableMemory;
			initialSize=(prealloc || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(mult*WAYS))) : initialSizeDefault);
			if(initialSize!=initialSizeDefault){
				System.err.println("Initial size set to "+initialSize);
			}
		}
		
		/* Adjust I/O settings and filenames */
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		for(int i=0; i<in1.size(); i++){
			String s=in1.get(i);
			if(s!=null && s.contains("#") && !new File(s).exists()){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				in1.set(i, a+1+b);
				in2.add(a+2+b);
			}
		}
		
		for(int i=0; i<ine1.size(); i++){
			String s=ine1.get(i);
			if(s!=null && s.contains("#") && !new File(s).exists()){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				ine1.set(i, a+1+b);
				ine2.add(a+2+b);
			}
		}
		
		for(int i=0; i<oute1.size(); i++){
			String s=oute1.get(i);
			if(s!=null && s.contains("#")){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				oute1.set(i, a+1+b);
				oute2.add(a+2+b);
			}
		}
		
//		if(in2!=null){
//			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
//			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
//		}

		if(!Tools.testOutputFiles(overwrite, append, false, outKmers, outContigs, outHist)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testOutputFiles(overwrite, append, false, oute1, oute2)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		if(!Tools.testInputFiles(true, true, ine1, ine2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		assert(THREADS>0);
		
		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			printMemory();
			outstream.println();
		}
		
		final int tableType=(useForest ? AbstractKmerTable.FOREST1D : useTable ? AbstractKmerTable.TABLE : useArray ? AbstractKmerTable.ARRAY1D : 0);
		tables=AbstractKmerTable.preallocate(WAYS, tableType, initialSize, (!prealloc || preallocFraction<1));
		
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void process(Timer t){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, outKmers, outHist);
		
		/* Count kmers */
		process2(processingMode);
		
		if(THREADS>1 && outHist!=null && outKmers!=null){
			Timer tout=new Timer();
			tout.start();
			Thread a=new DumpKmersThread();
			Thread b=new MakeKhistThread();
			a.start();
			b.start();
			while(a.getState()!=Thread.State.TERMINATED){
				try {
					a.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			while(b.getState()!=Thread.State.TERMINATED){
				try {
					b.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			tout.stop();
			outstream.println("Write Time:                 \t"+tout);
		}else{
			if(outHist!=null){
				makeKhist(outHist, histColumns, histMax, histHeader, histZeros, true);
			}
			if(outKmers!=null){
				dumpKmersAsText(outKmers, k, minToDump, true);
			}
		}
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		
		if(showSpeed){
			double rpnano=readsIn/(double)(t.elapsed);
			double bpnano=basesIn/(double)(t.elapsed);

			//Format with k or m suffixes
			String rpstring=(readsIn<100000 ? ""+readsIn : readsIn<100000000 ? (readsIn/1000)+"k" : (readsIn/1000000)+"m");
			String bpstring=(basesIn<100000 ? ""+basesIn : basesIn<100000000 ? (basesIn/1000)+"k" : (basesIn/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Total Time:               \t"+t);
			outstream.println("\nReads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		if(showStats && outContigs!=null && processingMode==contigMode && FileFormat.isFasta(ReadWrite.rawExtension(outContigs))){
			outstream.println();
			jgi.AssemblyStats2.main(new String[] {"in="+outContigs});
		}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException("BBDuk terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(int mode){
		
		/* Start phase timer */
		Timer t=new Timer();
		t.start();

		if(DISPLAY_PROGRESS){
			outstream.println("Before Load:");
			printMemory();
			outstream.println();
		}
		
		prefilterArray=null;
//		outstream.println();
		if(prefilter){
			KmerCountAbstract.CANONICAL=true;
			
			long precells=-1;
			int cbits=1;
			if(onePass){
				while(filterMax>=(1<<cbits)){cbits*=2;}
			}else{
				while(filterMax+1>=(1<<cbits)){cbits*=2;}
			}
			byte minq=0;
			if(precells<1){
				long prebits=(filterMemory-10)*8;
				precells=prebits/cbits;
			}
			if(prehashes<1){prehashes=2;}
			
			if(onePass){
				prefilterArray=KmerCount7MTA.makeKca(null, null, null, k, cbits, 0, precells, prehashes, minq, true, ecc, maxReads, 1, 1, 1, 1, null);
			}else if(precells>100000){
				Timer ht=new Timer();
				ht.start();
				
				ArrayList<String> extra=null;
				prefilterArray=KmerCount7MTA.makeKca_als(in1, in2, extra, k, cbits, 0, precells, prehashes, minq, true, ecc, maxReads, 1, 1, 1, 1, null);
				assert(filterMax<prefilterArray.maxValue);
				outstream.println("Made prefilter:   \t"+prefilterArray.toShortString(prehashes));
				double uf=prefilterArray.usedFraction();
				if(uf>0.6){
					outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" : 
						uf>0.7 ? "fairly" : "somewhat")+" full.  Ideal load is under 60% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; or increase prefiltersize.");
				}
				ht.stop();

				if(DISPLAY_PROGRESS){
					outstream.println("Prefilter time:\t"+ht);
					outstream.println("After prefilter:");
					printMemory();
					outstream.println();
				}
			}else{prefilter=false;}
		}
		
		/* Fill tables with kmers */
		long added=loadKmers();
		
		/* Clear prefilter; no longer needed */
		prefilterArray=null;
		
		if(DISPLAY_PROGRESS){
			outstream.println("\nAfter loading:");
			printMemory();
			outstream.println();
		}
		
		t.stop();
		outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		outstream.println("Unique Kmers:               \t"+added);
		outstream.println("Load Time:                  \t"+t);
		
		
		t.start();
		
		final long oldBasesIn=basesIn, oldReadsIn=readsIn;
		readsIn=basesIn=0;
		
		if(mode==extendMode){
			extendReads();
			
			if(DISPLAY_PROGRESS){
				outstream.println("\nAfter building contigs:");
				printMemory();
				outstream.println();
			}
			
			t.stop();

			outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");
			outstream.println("Output:                     \t"+readsIn+" reads \t\t"+(basesIn+basesExtended)+" bases.");
			outstream.println("Bases extended:             \t"+basesExtended);
			outstream.println("Reads extended:             \t"+readsExtended+String.format(" \t(%.2f%%)", readsExtended*100.0/readsIn));
		}else{
			/* Build contigs */
			buildContigs(mode);
			
			if(DISPLAY_PROGRESS){
				outstream.println("\nAfter building contigs:");
				printMemory();
				outstream.println();
			}
			
			t.stop();
			
			if(readsIn>0){outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");}
			outstream.println("Bases generated:            \t"+basesBuilt);
			outstream.println("Contigs generated:          \t"+contigsBuilt);
			outstream.println("Longest contig:             \t"+longestContig);
			outstream.println("Contig-building time:       \t"+t);
		}
		
		readsIn+=oldReadsIn;
		basesIn+=oldBasesIn;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private long loadKmers(){
		long added=0;
		for(int i=0; i<in1.size(); i++){
			String a=in1.get(i);
			String b=in2.size()>i ? in2.get(i) : null;
			int idx=a.indexOf('#');
			if(idx>=0 && b==null){
				b=a.replaceFirst("#", "2");
				a=a.replaceFirst("#", "1");
			}
			added+=loadKmers(a, b);
		}
		return added;
	}
	
	/**
	 * Load reads into tables, using multiple LoadThread.
	 */
	private long loadKmers(String fname1, String fname2){
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff1, ff2);
			cris.start(); //4567
		}
		
		/* Create ProcessThreads */
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new LoadThread(cris));}
		for(LoadThread pt : alpt){pt.start();}
		
		long added=0;
		
		/* Wait for threads to die, and gather statistics */
		for(LoadThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			added+=pt.added;
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
		return added;
	}
	
	
	/**
	 * Build contigs.
	 */
	private void buildContigs(final int mode){
		
		if(mode==contigMode){
			allContigs=new ArrayList<Read>();
			allInserts=null;
			
			if(useOwnership){
				for(AbstractKmerTable akt : tables){
					akt.initializeOwnership();
				}
			}
			
		}else if(mode==insertMode){
			allContigs=null;
			allInserts=new LongList();
		}else if(mode==extendMode){
			throw new RuntimeException("extendMode: TODO");
		}else{
			throw new RuntimeException("Unknown mode "+mode);
		}
		
		/* Create read input stream */
		final ConcurrentReadInputStream[] crisa=makeCrisArray(ine1, ine2);
		
		/* Create ProcessThreads */
		ArrayList<BuildThread> alpt=new ArrayList<BuildThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new BuildThread(i, mode, crisa));}
		for(BuildThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(BuildThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			for(Read contig : pt.contigs){
				allContigs.add(contig);
				contigsBuilt++;
				basesBuilt+=contig.length();
				longestContig=Tools.max(longestContig, contig.length());
			}
			if(allInserts!=null){
				allInserts.add(pt.insertSizes);
			}
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
		}
		
		/* Shut down I/O streams; capture error status */
		for(ConcurrentReadInputStream cris : crisa){
			errorState|=ReadWrite.closeStreams(cris);
		}
		
		if(outInsert!=null){
			FileFormat ff=FileFormat.testOutput(outInsert, FileFormat.TEXT, 0, 0, true, overwrite, append, false);
			TextStreamWriter tsw=new TextStreamWriter(ff);
			tsw.start();
			for(int i=0; i<allInserts.size; i++){
				long count=allInserts.get(i);
				if(count>0 || histZeros){
					tsw.print(i+"\t"+count+"\n");
				}
			}
			errorState|=tsw.poisonAndWait();
		}
		
		if(outContigs!=null){
			FileFormat ff=FileFormat.testOutput(outContigs, FileFormat.FA, 0, 0, true, overwrite, append, false);
//			ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, null, null, 4, null, false);
//			ros.start();
			ByteStreamWriter bsw=new ByteStreamWriter(ff);
			bsw.start();
			if(allContigs!=null){
//				Collections.sort(allContigs, ReadComparatorID.comparator);
				Collections.sort(allContigs, ReadLengthComparator.comparator);
				for(int i=0; i<allContigs.size(); i++){
					Read r=allContigs.get(i);
					bsw.println(r);
				}
			}
			errorState|=bsw.poisonAndWait();
		}
	}
	
	
	/**
	 * Extend reads.
	 */
	private void extendReads(){

		/* Create read input stream */
		final ConcurrentReadInputStream[] crisa=makeCrisArray(ine1, ine2);

		/* Create read input stream */
		final ConcurrentReadOutputStream[] rosa=makeCrosArray(oute1, oute2);
		
		/* Create ProcessThreads */
		ArrayList<ExtendThread> alpt=new ArrayList<ExtendThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ExtendThread(i, crisa, rosa));}
		for(ExtendThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ExtendThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
			readsExtended+=pt.readsExtendedT;
			basesExtended+=pt.basesExtendedT;
		}
		
		/* Shut down I/O streams; capture error status */
		for(ConcurrentReadInputStream cris : crisa){
			errorState|=ReadWrite.closeStreams(cris);
		}
		/* Shut down I/O streams; capture error status */
		if(rosa!=null){
			for(ConcurrentReadOutputStream ros : rosa){
				errorState|=ReadWrite.closeStream(ros);
			}
		}
	}
	
	private ConcurrentReadInputStream[] makeCrisArray(ArrayList<String> list1, ArrayList<String> list2){
		final ConcurrentReadInputStream[] array;

		array=new ConcurrentReadInputStream[list1.size()];
		for(int i=0; i<list1.size(); i++){
			String a=list1.get(i);
			String b=(list2.size()>i ? list2.get(i): null);
			if(verbose){System.err.println("Creating cris for "+a);}

			final ConcurrentReadInputStream cris;
			{
				FileFormat ff1=FileFormat.testInput(a, FileFormat.FASTA, null, true, true);
				FileFormat ff2=(b==null ? null : FileFormat.testInput(b, FileFormat.FASTA, null, true, true));
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff1, ff2);
			}
			array[i]=cris;
		}
		return array;
	}
	
	private ConcurrentReadOutputStream[] makeCrosArray(ArrayList<String> list1, ArrayList<String> list2){
		final ConcurrentReadOutputStream[] array;

		array=new ConcurrentReadOutputStream[list1.size()];
		for(int i=0; i<list1.size(); i++){
			String a=list1.get(i);
			String b=(list2.size()>i ? list2.get(i): null);
			if(verbose){System.err.println("Creating cris for "+a);}

			final ConcurrentReadOutputStream cris;
			{
				final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
				FileFormat ff1=FileFormat.testOutput(a, FileFormat.FASTQ, null, true, overwrite, append, ordered);
				FileFormat ff2=(b==null ? null : FileFormat.testOutput(b, FileFormat.FASTQ, null, true, overwrite, append, ordered));
				cris=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			}
			array[i]=cris;
		}
		return array;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads kmers. 
	 */
	private class LoadThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 */
		public LoadThread(ConcurrentReadInputStream cris_){
			cris=cris_;
			table=new HashBuffer(tables, buflen);
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					readsInT++;
					basesInT+=r1.length();
					if(r2!=null){
						readsInT++;
						basesInT+=r2.length();
					}
					
					if(ecc && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){BBMerge.findOverlapStrict(r1, r2, true);}

					if(r1!=null){
						if(r1.discarded()){
							lowqBasesT+=r1.length();
							lowqReadsT++;
						}else{
							long temp=addKmersToTable(r1);
							added+=temp;
							if(verbose){System.err.println("Added "+temp);}
						}
					}
					if(r2!=null){
						if(r2.discarded()){
							lowqBasesT+=r2.length();
							lowqReadsT++;
						}else{
							long temp=addKmersToTable(r2);
							added+=temp;
							if(verbose){System.err.println("Added "+temp);}
						}
					}
				}
				
				//Fetch a new read list
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			added+=table.flush();
		}
		
		
		private final int addKmersToTable(final Read r){
			if(onePass){return addKmersToTable_onePass(r);}
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final long key=toValue(kmer, rkmer);
					if(!prefilter || prefilterArray.read(key)>filterMax){
						int temp=table.incrementAndReturnNumCreated(key);
						created+=temp;
						if(verbose){System.err.println("Added "+temp);}
					}
				}
			}
			return created;
		}
		
		
		private final int addKmersToTable_onePass(final Read r){
			assert(prefilter);
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=~((-1L)<<shift);
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final long key=toValue(kmer, rkmer);
					int count=prefilterArray.incrementAndReturnUnincremented(key, 1);
					if(count>=filterMax){
						int temp=table.incrementAndReturnNumCreated(key);
						created+=temp;
						if(verbose){System.err.println("Added "+temp);}
					}
				}
			}
			return created;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream cris;
		
		private final HashBuffer table;
		
		public long added=0;
		
		private long readsInT=0;
		private long basesInT=0;
		private long lowqReadsT=0;
		private long lowqBasesT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds contigs. 
	 */
	private class BuildThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 */
		public BuildThread(int id_, int mode_, ConcurrentReadInputStream[] crisa_){
			id=id_;
			crisa=crisa_;
			mode=mode_;
		}
		
		@Override
		public void run(){
			if(crisa==null || crisa.length==0){
				//Build from kmers
				for(int tnum=id; tnum<tables.length; tnum+=THREADS){
//					System.err.println("id="+id+" processing table "+tnum);
					final HashArray1D table=(HashArray1D)tables[tnum];
					final int max=table.arrayLength();
					for(int cell=0; cell<max; cell++){
						int x=processCell(table, cell);
					}
				}
			}else{
				//Extend reads
				for(ConcurrentReadInputStream cris : crisa){
					synchronized(crisa){
						if(!cris.started()){
							cris.start();
						}
					}
					run(cris);
				}
			}
		}
		
		private int processCell(HashArray1D table, int cell){
			int count=table.readCellValue(cell);
			if(count<minCount){return 0;}
			
			long key=table.getKmer(cell);

			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+key+"\t"+AminoAcid.kmerToString(key, k));}
			if(useOwnership){
				int owner=table.getCellOwner(cell);
				if(verbose){outstream.println("Owner is initially "+owner);}
				if(owner>-1){return 0;}
				owner=table.setOwner(key, id, cell);
				if(verbose){outstream.println("Owner is now "+owner);}
				if(owner!=id){return 0;}
			}
			builderT.setLength(0);
			builderT.appendKmer(key, k);
			if(verbose){outstream.println("Filled builder: "+builderT);}
			
			byte[] contig=makeContig(builderT, id, true);
			if(contig!=null){
				if(verbose){System.err.println("Added "+contig.length);}
				final long num=contigNum.incrementAndGet();
				Read r=new Read(contig, -1, -1, -1, "*", null, num, 0);
				float gc=r.gc();
				float coverage=calcCoverage(contig, contig.length);
				r.id="contig_"+num+",length="+contig.length+",cov="+String.format("%.1f", coverage)+",gc="+String.format("%.3f", gc);
				contigs.add(r);
				return contig.length;
			}else{
				if(verbose){System.err.println("Created null contig.");}
			}
			return 0;
		}
		
		private void run(ConcurrentReadInputStream cris){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					processReadPair(r1, r2);
				}
				
				//Fetch a new read list
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}
		
		private void processReadPair(Read r1, Read r2){
			if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
			
			readsInT++;
			basesInT+=r1.length();
			if(r2!=null){
				readsInT++;
				basesInT+=r2.length();
			}
			
			if(mode==insertMode){
				int x=BBMerge.findOverlapStrict(r1, r2, false);
				if(x<1){
					x=findInsertSize(r1, r2, rightCounts);
				}
				insertSizes.increment(Tools.max(x, 0));
				return;
			}
			
			if(ecc && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){BBMerge.findOverlapStrict(r1, r2, true);}

			if(r1!=null){
				if(r1.discarded()){
					lowqBasesT+=r1.length();
					lowqReadsT++;
				}else{
					byte[] contig=makeContig(r1.bases, builderT, r1.numericID);
					if(contig!=null){
						if(verbose){System.err.println("Added "+contig.length);}
						final long num=contigNum.incrementAndGet();
						Read temp=new Read(contig, -1, -1, -1, "contig_"+num+"_length_"+contig.length, null, num, 0);
						contigs.add(temp);
					}
				}
			}
			if(r2!=null){
				if(r2.discarded()){
					lowqBasesT+=r2.length();
					lowqReadsT++;
				}else{
					byte[] contig=makeContig(r2.bases, builderT, r1.numericID);
					if(contig!=null){
						if(verbose){System.err.println("Added "+contig.length);}
						final long num=contigNum.incrementAndGet();
						Read temp=new Read(contig, -1, -1, -1, "contig_"+num+"_length_"+contig.length, null, num, 0);
						contigs.add(temp);
					}
				}
			}
		}
		
		/** From kmers */
		private byte[] makeContig(final ByteBuilder bb, long rid, boolean alreadyClaimed){
			final int initialLength=bb.length();
			if(initialLength<k){return null;}
			
			boolean success=alreadyClaimed ? true : (useOwnership ? claim(bb, id, rid) : true);
			if(verbose){System.err.println("Thread "+id+" checking owner after setting: "+findOwner(bb, id));}
			if(!success){
				release(bb, id);
				return null;
			}
			if(verbose  /*|| true*/){System.err.println("Thread "+id+" building contig; initial length "+bb.length());}
			if(verbose){System.err.println("Extending to right.");}
			success=extendToRight(bb, leftCounts, rightCounts, id);
			if(!success){
				release(bb, id);
				return null;
			}
			bb.reverseComplementInPlace();
			if(verbose  /*|| true*/){System.err.println("Extending rcomp to right; current length "+bb.length());}
//			verbose=true;
			success=extendToRight(bb, leftCounts, rightCounts, id);
			if(!success){
				release(bb, id);
				return null;
			}
			if(verbose  /*|| true*/){System.err.println("Final length for thread "+id+": "+bb.length());}
			if(bb.length()>initialLength+minExtension && bb.length()>=minContigLen){
//				if(useOwnership && THREADS==1){assert(claim(bases, bases.length, id, rid));}
				success=(useOwnership ? doubleClaim(bb, id, rid) : true);
				if(verbose  /*|| true*/){System.err.println("Success for thread "+id+": "+success);}
				if(success){
					bb.reverseComplementInPlace();
					return bb.toBytes();
				}else{
//					assert(false) : bb.length()+", "+id;
					release(bb.array, bb.length(), id);
					return null;
				}
			}else{
				success=false;
			}
			if(verbose  /*|| true*/){System.err.println("Contig was too short for "+id+": "+bb.length());}
			return null;
		}
		
		/** From a seed */
		private byte[] makeContig(final byte[] bases, final ByteBuilder bb, long rid){
			if(bases==null || bases.length<k){return null;}
//			if(verbose  /*|| true*/){System.err.println("Thread "+id+" checking owner: "+findOwner(bases, bases.length, id));}
			int owner=useOwnership ? findOwner(bases, bases.length, id) : -1;
			if(owner>=id){return null;}
			boolean success=(useOwnership ? claim(bases, bases.length, id, rid) : true);
			if(verbose  /*|| true*/){System.err.println("Thread "+id+" checking owner after setting: "+findOwner(bases, bases.length, id));}
			if(!success){
				release(bases, bases.length, id);
				return null;
			}
			if(verbose  /*|| true*/){System.err.println("Thread "+id+" building contig; initial length "+bases.length);}
			bb.setLength(0);
			bb.append(bases);
			if(verbose){System.err.println("Extending to right.");}
			success=extendToRight(bb, leftCounts, rightCounts, id);
			if(!success){
				release(bb.array, bb.length(), id);
				return null;
			}
			bb.reverseComplementInPlace();
			if(verbose  /*|| true*/){System.err.println("Extending rcomp to right; current length "+bb.length());}
//			verbose=true;
			success=extendToRight(bb, leftCounts, rightCounts, id);
			if(!success){
				release(bb.array, bb.length(), id);
				return null;
			}
			if(verbose  /*|| true*/){System.err.println("Final length for thread "+id+": "+bb.length());}
			if(bb.length()>bases.length+minExtension && bb.length()>=minContigLen){
//				if(useOwnership && THREADS==1){assert(claim(bases, bases.length, id, rid));}
				success=(useOwnership ? doubleClaim(bb, id, rid) : true);
				if(verbose  /*|| true*/){System.err.println("Success for thread "+id+": "+success);}
				if(success){
					bb.reverseComplementInPlace();
					return bb.toBytes();
				}else{
//					assert(false) : bb.length()+", "+id;
					release(bb.array, bb.length(), id);
					return null;
				}
			}else{
				success=false;
			}
			if(verbose  /*|| true*/){System.err.println("Contig was too short for "+id+": "+bb.length());}
			return null;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream[] crisa;
		
		private final int mode;

		private final int[] leftCounts=new int[4];
		private final int[] rightCounts=new int[4];
		private final ByteBuilder builderT=new ByteBuilder();
		
		private final LongList insertSizes=new LongList();
		
		ArrayList<Read> contigs=new ArrayList<Read>();
		
		private long readsInT=0;
		private long basesInT=0;
		private long lowqReadsT=0;
		private long lowqBasesT=0;
		private final int id;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Extends reads. 
	 */
	private class ExtendThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 */
		public ExtendThread(int id_, ConcurrentReadInputStream[] crisa_, ConcurrentReadOutputStream[] rosa_){
			id=id_;
			crisa=crisa_;
			rosa=rosa_;
			leftCounts=extendThroughLeftJunctions ? null : new int[4];
		}
		
		@Override
		public void run(){
			for(int i=0; i<crisa.length; i++){
				ConcurrentReadInputStream cris=crisa[i];
				ConcurrentReadOutputStream ros=(rosa!=null && rosa.length>i ? rosa[i] : null);
				synchronized(crisa){
					if(!cris.started()){
						cris.start();
					}
				}
				if(ros!=null){
					synchronized(rosa){
						if(!ros.started()){
							ros.start();
						}
					}
				}
				run(cris, ros);
			}
		}
		
		private void run(ConcurrentReadInputStream cris, ConcurrentReadOutputStream ros){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				final ArrayList<Read> listOut=new ArrayList<Read>(reads.size());
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					processReadPair(r1, r2);
					listOut.add(r1);
				}
				if(ros!=null){ros.add(listOut, ln.id);}
				
				//Fetch a new read list
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}
		
		private void processReadPair(Read r1, Read r2){
			if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
			
			readsInT++;
			basesInT+=r1.length();
			if(r2!=null){
				readsInT++;
				basesInT+=r2.length();
			}
			
			if(ecc && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){BBMerge.findOverlapStrict(r1, r2, true);}

			if(r1!=null){
				if(r1.discarded()){
					lowqBasesT+=r1.length();
					lowqReadsT++;
				}else{
					int extension=0;
					if(extendRight>0){
						extension+=extendRead(r1, builderT, leftCounts, rightCounts, extendRight);
					}
					if(extendLeft>0){
						r1.reverseComplement();
						extension+=extendRead(r1, builderT, leftCounts, rightCounts, extendLeft);
						r1.reverseComplement();
					}
					basesExtendedT+=extension;
					readsExtendedT+=(extension>0 ? 1 : 0);
				}
			}
			if(r2!=null){
				if(r2.discarded()){
					lowqBasesT+=r2.length();
					lowqReadsT++;
				}else{
					int extension=0;
					if(extendRight>0){
						extension+=extendRead(r2, builderT, leftCounts, rightCounts, extendRight);
					}
					if(extendLeft>0){
						r2.reverseComplement();
						extension+=extendRead(r2, builderT, leftCounts, rightCounts, extendLeft);
						r2.reverseComplement();
					}
					basesExtendedT+=extension;
					readsExtendedT+=(extension>0 ? 1 : 0);
				}
			}
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream[] crisa;
		private final ConcurrentReadOutputStream[] rosa;

		private final int[] leftCounts;
		private final int[] rightCounts=new int[4];
		private final ByteBuilder builderT=new ByteBuilder();
		
		private long readsInT=0;
		private long basesInT=0;
		private long lowqReadsT=0;
		private long lowqBasesT=0;
		private long readsExtendedT=0;
		private long basesExtendedT=0;
		private final int id;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Extension Methods      ----------------*/
	/*--------------------------------------------------------------*/


	public int findInsertSize(Read r1, Read r2, int[] rightCounts){
		final long kmer1=rightmostKmer(r1.bases, r1.length());
		final long kmer2=rightmostKmer(r2.bases, r2.length());
		if(kmer1<0 || kmer2<0){return -1;}
		final long rkmer1=AminoAcid.reverseComplementBinaryFast(kmer1, k);
		final long rkmer2=AminoAcid.reverseComplementBinaryFast(kmer2, k);
		final int x=measureInsert(kmer1, rkmer1, kmer2, rkmer2, 24000, rightCounts);
		if(x<0){return -1;}
		return r1.length()+r2.length()+x-k;//TODO: May be off by 1.
	}

	public int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance){
		final int initialLen=r.length();
		if(initialLen<k){return 0;}
		bb.setLength(0);
		bb.append(r.bases);
		final int extension=extendToRight2(bb, leftCounts, rightCounts, distance, true);
		if(extension>0){
			r.bases=bb.toBytes();
			if(r.quality!=null){
				r.quality=Arrays.copyOf(r.quality, r.bases.length);
				for(int i=initialLen; i<r.quality.length; i++){
					r.quality[i]=20;
				}
			}
		}
		assert(extension==r.length()-initialLen);
		return extension;
	}
	
	
	/** Returns distance between the two kmers, or -1 */
	public int measureInsert(final long kmer1, final long rkmer1, final long kmer2, final long rkmer2, final int maxlen, final int[] rightCounts){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		final long key2=toValue(kmer2, rkmer2);
		long kmer=kmer1;
		long rkmer=rkmer1;
		int len=0;
		
		{
			int way=(int)(key2%WAYS);
			int count=tables[way].getValue(key2);
			if(count<minCount){return -1;}
		}
		
		long key=toValue(kmer, rkmer);
		int way=(int)(key%WAYS);
		int count=tables[way].getValue(key);
		if(count<minCount){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return -1;
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
//		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
//		int rightSecond=rightCounts[rightSecondPos];
		
		if(rightMax<minCount){return -1;}
//		if(isJunction(rightMax, rightSecond)){return -1;}
		
		while(key!=key2 && len<maxlen){
			
			//Generate the new kmer
//			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			
			key=toValue(kmer, rkmer);
			way=(int)(key%WAYS);
			
			assert(tables[way].getValue(key)==rightMax);
			count=rightMax;
			
			assert(count>=minCount) : count;
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
			rightMax=rightCounts[rightMaxPos];
//			rightSecondPos=Tools.secondHighestPosition(rightCounts);
//			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+AbstractKmerTable.toText(kmer, k)+", "+AbstractKmerTable.toText(rkmer, k));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
//				outstream.println("rightSecondPos="+rightSecondPos);
//				outstream.println("rightSecond="+rightSecond);
			}
			
			if(rightMax<minCount){
				if(verbose){outstream.println("Breaking because highest right was too low:"+rightMax);}
				break;
			}

//			if(isJunction(rightMax, rightSecond)){return -1;}
			
			len++;
		}
		return len>=maxlen ? -1 : len;
	}
	

	
	/**
	 * Extend these bases into a contig.
	 * Stops at both left and right junctions.
	 * Claims ownership.
	 */
	public boolean extendToRight(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int id){
		if(bb.length()<k){return false;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			final int bblen=bb.length();
			final byte[] bases=bb.array;
			for(int i=bblen-k; i<bblen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}
		
		if(len<k){return false;}
		else{assert(len==k);}
		
		/* Now the trailing kmer has been initialized. */
		
		long key=toValue(kmer, rkmer);
		int way=(int)(key%WAYS);
		int count=tables[way].getValue(key);
		if(count<minCount){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return false;
		}
		
		int owner=(useOwnership ? tables[way].getOwner(key) : id);
		if(verbose){outstream.println("Owner: "+owner);}
		if(owner>id){return false;}
		
		int leftMaxPos=0;
		int leftMax=minCount;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+AbstractKmerTable.toText(kmer, k)+", "+AbstractKmerTable.toText(rkmer, k));
			outstream.println("Counts: "+count+", "+(leftCounts==null ? "null" : Arrays.toString(leftCounts))+", "+Arrays.toString(rightCounts));
			outstream.println("leftMaxPos="+leftMaxPos);
			outstream.println("leftMax="+leftMax);
			outstream.println("leftSecondPos="+leftSecondPos);
			outstream.println("leftSecond="+leftSecond);
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}
		
		if(rightMax<minCount){return false;}
		if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){return false;}
		
		if(useOwnership){
			owner=tables[way].setOwner(key, id);
			if(verbose){outstream.println("A. Owner is now "+id+" for key "+key);}
			if(owner!=id){
				if(verbose){outstream.println("Returning early because owner was "+owner+" for thread "+id+".");}
				return false;
			}
		}
		
		final int maxLen=Tools.min((extendRight<0 ? maxContigLen : bb.length()+extendRight), maxContigLen);
		
		while(owner==id && bb.length()<maxLen){
			
			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			
			key=toValue(kmer, rkmer);
			way=(int)(key%WAYS);
			
			assert(tables[way].getValue(key)==rightMax);
			count=rightMax;
			
			assert(count>=minCount) : count;

			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+AbstractKmerTable.toText(kmer, k)+", "+AbstractKmerTable.toText(rkmer, k));
				outstream.println("Counts: "+count+", "+(leftCounts==null ? "null" : Arrays.toString(leftCounts))+", "+Arrays.toString(rightCounts));
				outstream.println("leftMaxPos="+leftMaxPos);
				outstream.println("leftMax="+leftMax);
				outstream.println("leftSecondPos="+leftSecondPos);
				outstream.println("leftSecond="+leftSecond);
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			
			if(rightMax<minCount){
				if(verbose){outstream.println("Breaking because highest right was too low:"+rightMax);}
				break;
			}
			
			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){break;}
			
			if(useOwnership){
				owner=tables[way].getOwner(key);
				if(verbose){outstream.println("Owner is initially "+id+" for key "+key);}
				if(owner==id){
					if(verbose  /*|| true*/){
//						outstream.println(new String(bb.array, bb.length()-31, 31));
						outstream.println(bb);
						outstream.println(AbstractKmerTable.toText(kmer, 31));
						outstream.println(AbstractKmerTable.toText(rkmer, 31));
						outstream.println("Breaking because owner was "+owner+" for thread "+id+".");
					}
					break;
				}
				owner=tables[way].setOwner(key, id);
				if(verbose){outstream.println("B. Owner is now "+id+" for key "+key);}
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		if(verbose  /*|| true*/){outstream.println("Returning because owner was "+owner+" for thread "+id+".");}
		return id==owner;
	}
	

	
	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 */
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase){
		final int initialLength=bb.length();
		if(initialLength<k){return 0;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0;
		long rkmer=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			int len=0;
			final byte[] bases=bb.array;
			for(int i=initialLength-k; i<initialLength; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
			if(len<k){return 0;}
			else{assert(len==k);}
		}
		
		/* Now the trailing kmer has been initialized. */
		
		long key=toValue(kmer, rkmer);
		int way=(int)(key%WAYS);
		int count=tables[way].getValue(key);
		if(count<minCount){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return 0;
		}
		
		int leftMaxPos=0;
		int leftMax=minCount;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+AbstractKmerTable.toText(kmer, k)+", "+AbstractKmerTable.toText(rkmer, k));
			outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}
		
		if(rightMax<minCount){return 0;}
		if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){return 0;}
		
		final int maxLen=Tools.min(bb.length()+distance, maxContigLen);
		
		while(bb.length()<maxLen){
			
			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			
			key=toValue(kmer, rkmer);
			way=(int)(key%WAYS);
			
			assert(tables[way].getValue(key)==rightMax);
			count=rightMax;
			
			assert(count>=minCount) : count;
			
			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+AbstractKmerTable.toText(kmer, k)+", "+AbstractKmerTable.toText(rkmer, k));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			
			if(rightMax<minCount){
				if(verbose){outstream.println("Breaking because highest right was too low:"+rightMax);}
				break;
			}
			
			if(includeJunctionBase && kmer>rkmer){
				bb.append(b);
				if(verbose){outstream.println("Added base "+(char)b);}
			}

			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){break;}

			if(!includeJunctionBase || kmer<=rkmer){
				bb.append(b);
				if(verbose){outstream.println("Added base "+(char)b);}
			}
		}
		return bb.length()-initialLength;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private long rightmostKmer(final byte[] bases, final int blen){
		if(blen<k){return -1;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0;
		long rkmer=0;
		int len=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			for(int i=blen-k; i<blen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}
		
		if(len<k){return -1;}
		else{assert(len==k);}
		return kmer;
	}
	
	
	private boolean isJunction(int rightMax, int rightSecond, int leftMax, int leftSecond){
		if(isJunction(rightMax, rightSecond)){return true;}
		return isJunction(leftMax, leftSecond);
	}
	
	private boolean isJunction(int max, int second){
		if(second*branchMult1<max || (second<=branchLowerConst && max>=Tools.max(minCount, second*branchMult2))){
			return false;
		}
		if(verbose){outstream.println("Breaking because second-highest was too high:" +
				(second*branchMult1<max)+", "+(second<=branchLowerConst)+", "+(max>=Tools.max(minCount, second*branchMult2)));}
		return true;
	}
	
	private boolean doubleClaim(final ByteBuilder bb, final int id, final long rid){
		return doubleClaim(bb.array, bb.length(), id, rid);
	}
	
	/** Ensures there can be only one owner. */
	private boolean doubleClaim(final byte[] bases, final int blength, final int id, final long rid){
		boolean success=claim(bases, blength, id, rid);
		if(verbose){outstream.println("success1="+success+", id="+id+", blength="+blength);}
		if(!success){return false;}
		success=claim(bases, blength, id+CLAIM_OFFSET, rid);
		if(verbose){outstream.println("success2="+success+", id="+id+", blength="+blength);}
		return success;
	}
	
	private boolean claim(final ByteBuilder bb, final int id, final long rid){
		return claim(bb.array, bb.length(), id, rid);
	}
	
	private float calcCoverage(final byte[] bases, final int blength){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0, rkmer=0;
		int len=0;
		long sum=0, max=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(len>=k){
				int count=getCount(kmer, rkmer);
				sum+=count;
				max=Tools.max(count, max);
			}
		}
		return sum==0 ? 0 : sum/(float)blength;
	}
	
	private boolean claim(final byte[] bases, final int blength, final int id, final long rid){
		if(verbose){outstream.println("Thread "+id+" claim start.");}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0, rkmer=0;
		int len=0;
		boolean success=true;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength && success; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=k){
				success=claim(kmer, rkmer, id, rid, i);
			}
		}
		return success;
	}
	
	private boolean claim(final long kmer, final long rkmer, final int id, final long rid, final int pos){
		final long key=toValue(kmer, rkmer);
		final int way=(int)(key%WAYS);
		final AbstractKmerTable table=tables[way];
		final int count=table.getValue(key);
		assert(count==-1 || count>0) : count;
//		if(verbose  /*|| true*/){outstream.println("Count="+count+".");}
		if(count<0){return true;}
		assert(count>0) : count;
		final int owner=table.setOwner(key, id);
		if(verbose){outstream.println("owner="+owner+".");}
//		assert(owner==id) : id+", "+owner+", "+rid+", "+pos;
		return owner==id;
	}
	
	private void release(ByteBuilder bb, final int id){
		release(bb.array, bb.length(), id);
	}
	
	private void release(final byte[] bases, final int blength, final int id){
		if(verbose  /*|| true*/){outstream.println("*Thread "+id+" release start.");}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0, rkmer=0;
		int len=0;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=k){
				release(kmer, rkmer, id);
			}
		}
	}
	
	private boolean release(final long kmer, final long rkmer, final int id){
		final long key=toValue(kmer, rkmer);
		final int way=(int)(key%WAYS);
		final AbstractKmerTable table=tables[way];
		final int count=table.getValue(key);
//		if(verbose  /*|| true*/){outstream.println("Count="+count+".");}
		if(count<1){return true;}
		return table.clearOwner(key, id);
	}
	
	private int findOwner(ByteBuilder bb, final int id){
		return findOwner(bb.array, bb.length(), id);
	}
	
	private int findOwner(final byte[] bases, final int blength, final int id){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0, rkmer=0;
		int len=0;
		int maxOwner=-1;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=k){
				int owner=findOwner(kmer, rkmer);
				maxOwner=Tools.max(owner, maxOwner);
				if(maxOwner>id){break;}
			}
		}
		return maxOwner;
	}
	
	private int findOwner(final long kmer, final long rkmer){
		final long key=toValue(kmer, rkmer);
		final int way=(int)(key%WAYS);
		final AbstractKmerTable table=tables[way];
		final int count=table.getValue(key);
		if(count<0){return -1;}
		final int owner=table.getOwner(key);
		return owner;
	}

	private int getCount(long kmer, long rkmer){
		long key=toValue(kmer, rkmer);
		int way=(int)(key%WAYS);
		return tables[way].getValue(key);
	}
	
	private int fillRightCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==AminoAcid.reverseComplementBinary(rkmer, k));
		if(verbose){outstream.println("fillRightCounts:   "+AbstractKmerTable.toText(kmer, k)+",   "+AbstractKmerTable.toText(rkmer, k));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;
		
		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
			if(verbose){outstream.println("kmer:               "+AbstractKmerTable.toText(kmer2, k)+", "+AbstractKmerTable.toText(rkmer2, k));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==AminoAcid.reverseComplementBinary(rkmer2, k));
			long key=toValue(kmer2, rkmer2);
			int way=(int)(key%WAYS);
			int count=tables[way].getValue(key);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}
	
	private int fillLeftCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==AminoAcid.reverseComplementBinary(rkmer, k));
		if(verbose){outstream.println("fillLeftCounts:    "+AbstractKmerTable.toText(kmer, k)+",   "+AbstractKmerTable.toText(rkmer, k));}
		rkmer=(rkmer<<2)&mask;
		kmer=(kmer>>>2);
		int max=-1, maxPos=0;
//		assert(false) : shift2+", "+k;
		for(int i=0; i<=3; i++){
			long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
			long kmer2=kmer|(((long)i)<<shift2);
			if(verbose){outstream.println("kmer:             "+AbstractKmerTable.toText(kmer2, k)+",     "+AbstractKmerTable.toText(rkmer2, k));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==AminoAcid.reverseComplementBinary(rkmer2, k)) : "\n"+"kmer:      \t"+
				AbstractKmerTable.toText(AminoAcid.reverseComplementBinary(rkmer2, k), k)+", "+AbstractKmerTable.toText(AminoAcid.reverseComplementBinary(kmer2, k), k);
			long key=toValue(rkmer2, kmer2);
			int way=(int)(key%WAYS);
			int count=tables[way].getValue(key);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}
	
	public boolean dumpKmersAsText(String fname, int k, int minToDump, boolean printTime){
		if(fname==null){return false;}
		Timer t=new Timer();
		t.start();
		
//		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, true);
//		tsw.start();
//		for(AbstractKmerTable set : tables){
//			set.dumpKmersAsText(tsw, k, minToDump);
//		}
//		tsw.poisonAndWait();
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();
		for(AbstractKmerTable set : tables){
			set.dumpKmersAsBytes(bsw, k, minToDump);
		}
		bsw.poisonAndWait();
		
		t.stop();
		if(printTime){outstream.println("Kmer Dump Time:             \t"+t);}
		return bsw.errorState;
	}
	
	public boolean makeKhist(String fname, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime){
		if(fname==null){return false;}
		Timer t=new Timer();
		t.start();
		
//		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, true);
//		tsw.start();
//		if(printHeader){
//			tsw.print("#Depth\t"+(cols==3 ? "RawCount\t" : "")+"Count\n");
//		}
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();
		if(printHeader){
			bsw.print("#Depth\t"+(cols==3 ? "RawCount\t" : "")+"Count\n");
		}
		
		CoverageArray3 ca=new CoverageArray3(-1, histMax+1);
		for(AbstractKmerTable set : tables){
			set.fillHistogram(ca, histMax);
		}
		
//		StringBuilder sb=new StringBuilder();
		for(int i=1; i<=ca.maxIndex; i++){
			int count=ca.get(i);
			if(printZeros || count>0){
//				sb.append(i);
//				sb.append('\t');
//				if(cols==3){
//					sb.append(i*(long)count);
//					sb.append('\t');
//				}
//				sb.append(count);
//				sb.append('\n');
//				bsw.print(sb.toString());
//				sb.setLength(0);

				bsw.print(i);
				bsw.print('\t');
				if(cols==3){
					bsw.print(i*(long)count);
					bsw.print('\t');
				}
				bsw.print(count);
				bsw.print('\n');
			}
		}
		bsw.poisonAndWait();
		t.stop();
		if(printTime){outstream.println("Histogram Write Time:       \t"+t);}
		return bsw.errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class DumpKmersThread extends Thread {
		
		DumpKmersThread(){}
		
		public void run(){
			dumpKmersAsText(outKmers, k, minToDump, false);
		}
		
	}
	
	private class MakeKhistThread extends Thread {
		
		MakeKhistThread(){}
		
		public void run(){
			makeKhist(outHist, histColumns, histMax, histHeader, histZeros, false);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print statistics about current memory use and availability */
	private static final void printMemory(){
		if(GC_BEFORE_PRINT_MEMORY){
			System.gc();
			System.gc();
		}
		Runtime rt=Runtime.getRuntime();
		long mmemory=rt.maxMemory()/1000000;
		long tmemory=rt.totalMemory()/1000000;
		long fmemory=rt.freeMemory()/1000000;
		long umemory=tmemory-fmemory;
		outstream.println("Memory: "+"max="+mmemory+"m, total="+tmemory+"m, "+"free="+fmemory+"m, used="+umemory+"m");
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	private final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<Read> allContigs;
	private LongList allInserts;
	private long contigsBuilt=0;
	private long basesBuilt=0;
	private long longestContig=0;
	
	private boolean extendThroughLeftJunctions=true;
	
	private int processingMode=contigMode;

	private int extendLeft=-1;
	private int extendRight=-1;
	
	/** Track kmer ownership */
	public final boolean useOwnership;
	
//	public int maxExtension=500000000;
	public int maxContigLen=1000000000;
	public int minExtension=2;
	public int minContigLen=100;

	private int minCount=3;
	private float branchMult1=60;
	private float branchMult2=8;
	private int branchLowerConst=3;
	
	public boolean showStats=true;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;

	/** Use a count-min prefilter for low-depth kmers */
	public boolean prefilter=false;
	/** Fill the prefilter at the same time as the main table */
	public boolean onePass=false;
	/** Number of hashes used by prefilter */
	public int prehashes=2;
	/** Fraction of memory used by prefilter */
	private double prefilterFraction=0.2;
	
	/** Initial size of data structures */
	private int initialSize=-1;
	private static final int initialSizeDefault=128000;
	/** Fraction of available memory preallocated to arrays */
	private double preallocFraction=1.0;
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in tables[Y] */
	private final AbstractKmerTable[] tables;
	
	private KCountArray prefilterArray=null;

	/** Input reads for kmers */
	private ArrayList<String> in1=new ArrayList<String>(), in2=new ArrayList<String>();
	/** Input reads to extend */
	private ArrayList<String> ine1=new ArrayList<String>(), ine2=new ArrayList<String>();
	/** Output extended reads */
	private ArrayList<String> oute1=new ArrayList<String>(), oute2=new ArrayList<String>();
	
	
	/** Contig output file */
	private String outContigs=null;
	/** Insert size histogram */
	private String outInsert=null;
	/** Kmer count output file */
	private String outKmers=null;
	/** Histogram output file */
	private String outHist=null;

	/** Histogram columns */
	private int histColumns=2;
	/** Histogram rows */
	private int histMax=100000;
	/** Histogram columns */
	private boolean histHeader=false;
	/** Histogram show rows with 0 count */
	private boolean histZeros=false;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	
	private int buflen=1000;
	
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	long lowqReads=0;
	long lowqBases=0;
	long basesExtended=0;
	long readsExtended=0;
	
	private static final int CLAIM_OFFSET=100000;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	private final long usableMemory;
	private final long filterMemory;
	private final long tableMemory;
	
	/** Number of tables (and threads, during loading) */ 
	private final boolean prealloc;
	
	/** Number of tables (and threads, during loading) */ 
	private final int WAYS;
	
	/** Filter kmers up to this level; don't store them in primary data structure */ 
	private final int filterMax;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Use HashForest data structure */
	private final boolean useForest;
	/** Use KmerTable data structure */
	private final boolean useTable;
	/** Use HashArray data structure (default) */
	private final boolean useArray;	
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	
	/** min kmer count to dump to text */
	private int minToDump=1;
	
	/** Quality-trim the left side */
	private final boolean ecc;
	
	/** True iff java was launched with the -ea' flag */
	private final boolean EA;
	
	private AtomicLong contigNum=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=false;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Force output reads to stay in input order */
	public static boolean ordered=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Verbose messages */
	public static boolean verbose=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;

	private static final int contigMode=0;
	private static final int extendMode=1;
	private static final int insertMode=2;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/

	
}
