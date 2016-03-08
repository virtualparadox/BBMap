package jgi;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import bloom.KCountArray;
import bloom.KmerCount7MTA;
import bloom.KmerCountAbstract;


import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;

import dna.AminoAcid;
import dna.Data;
import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;

import align2.ListNum;
import align2.LongList;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;

/**
 * @author Brian Bushnell
 * @date Aug 14, 2012
 *
 */
public class BBMerge {
	
	
	public static void main(String[] args){
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		BBMerge mr=new BBMerge(args);
		mr.process();
		Read.VALIDATE_IN_CONSTRUCTOR=true;
	}
	
	
	private static void printOptions(){
		System.err.println("Syntax:\n");
		System.err.println("java -ea -Xmx200m -cp <path> jgi.BBMerge in=<input file> out=<output file>");
		System.err.println("For more options, please run the shellscript with no arguments, or look at its contents.");
	}
	
	
	private static String[] preparse(String[] args){
		if(args==null){return new String[0];}
		int nulls=0;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);}
			if(a.equals("jni") || a.equals("usejni")){
				Shared.USE_JNI=Tools.parseBoolean(b);
			}else if(a.equals("vstrict") || a.equals("verystrict")){
				vstrict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("ustrict") || a.equals("ultrastrict")){
				ustrict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("xstrict") || a.equals("hstrict") || a.equals("hyperstrict") || a.equals("maxstrict")){
				xstrict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("strict")){
				strict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("loose")){
				loose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("vloose") || a.equals("veryloose")){
				vloose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("uloose") || a.equals("ultraloose")){
				uloose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("xloose") || a.equals("hloose") || a.equals("hyperloose") || a.equals("maxloose")){
				xloose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("fast")){
				fast=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("default")){
				if(Tools.parseBoolean(b)){
					xstrict=ustrict=vstrict=strict=loose=vloose=uloose=xloose=fast=false;
				}
				args[i]=null;
				nulls++;
			}
		}
		
		if(nulls==0){return args;}
		ArrayList<String> args2=new ArrayList<String>(args.length-nulls+5);
		if(strict || vstrict || ustrict || xstrict){
			strict=true;
			loose=vloose=uloose=xloose=false;
			
			args2.add("maxbad=4");
			args2.add("margin=3");
			args2.add("minqo=8");
			args2.add("qualiters=2");
			
			if(xstrict){
				args2.add("ratiomode=t");
				args2.add("normalmode=t");
				args2.add("requireratiomatch=t");

				args2.add("minentropy=56");
				args2.add("minoverlap=14");
				args2.add("minoverlap0=3");
				
				args2.add("maxratio=0.055");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.65");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.25");
			}else if(ustrict){
				args2.add("ratiomode=t");
				args2.add("normalmode=t");
				args2.add("requireratiomatch=t");

				args2.add("minentropy=56");
				args2.add("minoverlap=14");
				args2.add("minoverlap0=3");
				
				args2.add("maxratio=0.045");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.5");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.03");
			}else if(vstrict){
				if(true){//faster mode
					args2.add("ratiomode=t");
					args2.add("normalmode=f");

					args2.add("minentropy=52");
					args2.add("minoverlap=12");
					args2.add("minoverlap0=4");

					args2.add("maxratio=0.05");
					args2.add("ratiomargin=12");
					args2.add("ratiooffset=0.5");
					args2.add("ratiominoverlapreduction=4");
					args2.add("efilter=2");
					args2.add("pfilter=0.008");
				}else{//slower but more accurate rrm mode
					args2.add("ratiomode=t");
					args2.add("normalmode=t");
					args2.add("requireratiomatch=t");

					args2.add("minentropy=42");
					args2.add("minoverlap=12");
					args2.add("minoverlap0=5");

					args2.add("maxratio=0.06");
					args2.add("ratiomargin=8");
					args2.add("ratiooffset=0.5");
					args2.add("ratiominoverlapreduction=4");
					args2.add("efilter=3");
				}
			}else{
				args2.add("ratiomode=t");
				args2.add("normalmode=f");
				
				args2.add("minentropy=42");
				args2.add("minoverlap0=7");
				args2.add("minoverlap=11");
				
				args2.add("maxratio=0.075");
				args2.add("ratiomargin=7.5");
				args2.add("ratiooffset=0.55");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=4");
				args2.add("pfilter=0.0008");
			}
		}else if(loose || vloose || uloose || xloose){
			loose=true;
			strict=vstrict=ustrict=xstrict=false;
			args2.add("minoverlap=8");
			args2.add("minoverlap0=9");
			args2.add("qualiters=4");
			args2.add("mismatches=3");
			args2.add("margin=2");
			
			args2.add("ratiooffset=0.4");
			
			if(xloose){
				args2.add("owq=t");
				args2.add("ouq=t");
				args2.add("minentropy=22");
				args2.add("minoverlap=8");
				args2.add("minoverlap0=7");
				args2.add("maxratio=0.2");
				args2.add("mismatches=3");
				args2.add("ratiomargin=2");
				args2.add("normalmode=t");
				args2.add("pfilter=0.0000001");
				args2.add("efilter=8");
				args2.add("margin=2");
				args2.add("ratiominoverlapreduction=2");
			}else if(vloose || uloose){
				args2.add("owq=t");
				args2.add("ouq=t");
				if(uloose){
//					args2.add("maxratio=0.14");
//					args2.add("ratiomargin=2");
//					args2.add("normalmode=t");
//					args2.add("pfilter=0.0000001");
					
					
					args2.add("minoverlap=8");
					args2.add("minoverlap0=7");
					args2.add("mismatches=3");
					args2.add("margin=2");

					args2.add("ratiominoverlapreduction=2");
					args2.add("efilter=8");
					args2.add("maxratio=0.16");
					args2.add("ratiomargin=2.2");
					args2.add("pfilter=0.0000002");
					args2.add("minentropy=24");
				}else{
					args2.add("ratiominoverlapreduction=3");
					args2.add("maxratio=0.12");
					args2.add("ratiomargin=3");
					args2.add("pfilter=0.000004");
					args2.add("minentropy=28");
					args2.add("efilter=7.5");
					args2.add("ratiooffset=0.45");
				}
			}else{
				args2.add("maxratio=0.11");
				args2.add("ratiomargin=4.7");
				args2.add("ratiominoverlapreduction=2");
				args2.add("pfilter=0.00002");
				args2.add("efilter=8");
				args2.add("minentropy=30");
			}
		}else if(fast){
			args2.add("maxratio=0.08");
			args2.add("ratiomargin=2.5");
			args2.add("ratiominoverlapreduction=3");
			args2.add("pfilter=0.0002");
			args2.add("efilter=8");
			args2.add("minentropy=39");
			args2.add("mininsert0=50");
		}
		
		for(String s : args){
			if(s!=null){args2.add(s);}
		}
		return args2.toArray(new String[args2.size()]);
	}
	
	public BBMerge(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		System.err.println("BBMerge version "+version);
		
		args=preparse(args);
		Timer ttotal=new Timer();
		ttotal.start();
		
		in1=(args[0].indexOf('=')>0 ? null : args[0]);
		in2=(in1!=null && args.length>1 && args[1].indexOf('=')<0 ? args[1] : null);
		if(in2!=null && "null".equalsIgnoreCase(in2)){in2=null;}
		
		{
			if(in1!=null && !in1.contains(",") && !in1.startsWith("stdin.") && !in1.equals("stdin")){
				File f=new File(in1);
				if(!f.exists() || !f.isFile()){
					in1=null;
//					throw new RuntimeException(in1+" does not exist.");
				}
			}
			if(in2!=null && !in2.contains(",")){
				File f=new File(in2);
				if(!f.exists() || !f.isFile()){
					in2=null;
//					throw new RuntimeException(in2+" does not exist.");
				}else if(in1.equalsIgnoreCase(in2)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		ReadWrite.MAX_ZIP_THREADS=Shared.threads()-1;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		Shared.READ_BUFFER_LENGTH=Tools.max(Shared.READ_BUFFER_LENGTH, 400);
		
		boolean mm0set=false;
		KmerCountAbstract.minProb=0f;
		
		Parser parser=new Parser();
		parser.trimq=trimq;
		parser.minAvgQuality=minAvgQuality;
		parser.minReadLength=minReadLength;
		parser.maxReadLength=maxReadLength;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);}

			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("useratio") || a.equals("ratio") || a.equals("ratiomode")){
				useRatioMode=Tools.parseBoolean(b);
			}else if(a.equals("usenormalmode") || a.equals("normalmode")){
				useNormalMode=Tools.parseBoolean(b);
			}else if(a.equals("requireratiomatch") || a.equals("rrm")){
				requireRatioMatch=Tools.parseBoolean(b);
			}else if(a.equals("maxratio")){
				MAX_RATIO=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiomargin")){
				RATIO_MARGIN=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiooffset")){
				RATIO_OFFSET=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiominoverlapreduction")){
				MIN_OVERLAPPING_BASES_RATIO_REDUCTION=Integer.parseInt(b);
//				useRatioMode=true;
			}else if(a.equals("minentropy") || a.equals("entropy")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					minEntropyScore=Integer.parseInt(b);
				}else{
					useEntropy=Tools.parseBoolean(b);
				}
			}else if(a.equals("minoverlappingbases") || a.equals("minoverlapbases") || a.equals("minoverlap")){
				MIN_OVERLAPPING_BASES=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases0") || a.equals("minoverlapbases0") || a.equals("minoverlap0")){
				MIN_OVERLAPPING_BASES_0=Integer.parseInt(b);
			}else if(a.equals("minqo") || a.equals("minq")){
				MIN_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("maxq")){
				Read.MAX_MERGE_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("qualiters")){
				QUAL_ITERS=Tools.max(1, Integer.parseInt(b));
			}else if(a.equals("maxbadbases") || a.equals("maxbad") || a.equals("mismatches")){
				MAX_MISMATCHES=Integer.parseInt(b);
			}else if(a.equals("maxbadbases0") || a.equals("maxbad0") || a.equals("mismatches0")){
				MAX_MISMATCHES0=Integer.parseInt(b);
				mm0set=true;
			}else if(a.equals("margin")){
				MISMATCH_MARGIN=Integer.parseInt(b);
			}else if(a.equals("usemapping")){
				USE_MAPPING=Tools.parseBoolean(b);
			}else if(a.equals("bin")){
				bin=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Shared.setThreads(Integer.parseInt(b));
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.equals("outgood") || a.equals("outmerged") || a.equals("outm") || a.equals("out")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood1") || a.equals("outmerged1") || a.equals("outm1") || a.equals("out1")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood2") || a.equals("outmerged2") || a.equals("outm2") || a.equals("out2")){
				out2=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb") || a.equals("outu") || a.equals("outunmerged") || a.equals("outbad")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb1") || a.equals("outu1") || a.equals("outunmerged1") || a.equals("outbad1")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb2") || a.equals("outu2") || a.equals("outunmerged2") || a.equals("outbad2")){
				outb2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outinsert") || a.startsWith("outlength")){
				outinsert=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist") || a.equals("hist") || a.equals("histogram") || a.equals("ihist")){
				outhist=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist2") || a.equals("hist2")){
				outhist2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist3") || a.equals("hist3")){
				outhist3=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outa") || a.equals("outadapter")){
				outAdapter=b;
				findAdapterSequence=(outAdapter!=null);
			}else if(a.equals("outputfailed")){
				OUTPUT_FAILED=Tools.parseBoolean(b);
			}else if(a.equals("mix")){
				MIX_BAD_AND_GOOD=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				NONZERO_ONLY=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "verbose flag is static final; recompile to change it.";
//				verbose=Tools.parseBoolean(b);
			}else if(a.equals("join") || a.equals("merge")){
				join=Tools.parseBoolean(b);
				if(join){ecc=false;}
			}else if(a.equals("ecc") || a.equals("errorcorrect")){
				ecc=Tools.parseBoolean(b);
				if(ecc){join=false;}
			}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
				trimByOverlap=Tools.parseBoolean(b);
			}else if(a.equals("useoverlap") || a.equals("usebases") || a.equals("matebyoverlap") || a.equals("matebybases")){
				MATE_BY_OVERLAP=Tools.parseBoolean(b);
			}else if(a.startsWith("skipmated")){
				SKIP_MATED_READS=Tools.parseBoolean(b);
			}else if(a.equals("lowercase")){
				lowercaseAdapters=Tools.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("trimonfailure") || a.equals("tof")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					TRIM_ON_OVERLAP_FAILURE=Integer.parseInt(b);
				}else{
					TRIM_ON_OVERLAP_FAILURE=(Tools.parseBoolean(b) ? 1 : 0);
				}
			}else if(a.equals("overlapusingquality") || a.equals("ouq")){
				overlapUsingQuality=Tools.parseBoolean(b);
			}else if(a.equals("overlapwithoutquality") || a.equals("owoq") || a.equals("owuq") || a.equals("owq")){
				overlapWithoutQuality=Tools.parseBoolean(b);
			}else if(a.equals("maxExpectedErrors") || a.equals("mee") || a.equals("meefilter")){
				maxExpectedErrors=Float.parseFloat(b);
			}else if(a.equals("mi") || a.equals("minins") || a.equals("mininsert")){
				minInsert=Integer.parseInt(b);
			}else if(a.equals("mi0") || a.equals("mininsert0")){
				minInsert0=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				KmerCountAbstract.minProb=Float.parseFloat(b);
				assert(KmerCountAbstract.minProb<1) : "minprob must be less than 1.  At 1, even kmers with 100% probablity of being error-free will be discarded.";
			}else if(a.equals("minqf") || a.equals("minqfilter")){
				filterMinq=(byte)Integer.parseInt(b);
			}else if(a.equals("k")){
				filterK=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				filterHashes=Integer.parseInt(b);
			}else if(a.equals("efilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					boolean x=Tools.parseBoolean(b);
					if(!x){efilterRatio=0;}
				}else{
					efilterRatio=Float.parseFloat(b);
				}
				useEfilter=efilterRatio>0;
			}else if(a.equals("pfilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					boolean x=Tools.parseBoolean(b);
					if(!x){pfilterRatio=0;}
				}else{
					pfilterRatio=Float.parseFloat(b);
				}
			}else if(a.equals("efilteroffset")){
				efilterOffset=Float.parseFloat(b);
			}else if(a.equals("kfilter")){
				useKFilter=Tools.parseBoolean(b);
			}else if(a.equals("usequality")){
				useQuality=Tools.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Tools.parseBoolean(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("recalibrate") || a.equals("recalibratequality") || a.equals("recal")){
				recalibrateQuality=Tools.parseBoolean(b);
			}else if(a.equals("recalpairnum") || a.equals("recalibratepairnum")){
				CalcTrueQuality.USE_PAIRNUM=Tools.parseBoolean(b);
			}else if(a.equals("path")){
				Data.setPath(b);
			}else if(a.equals("iupacton") || a.equals("itn")){
				iupacToN=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		minInsert=Tools.max(minInsert, MIN_OVERLAPPING_BASES);
		if(minInsert0<1){
			minInsert0=(Tools.max((int)(minInsert*0.75), 5, MIN_OVERLAPPING_BASES_0));
			int cap=(loose ? 50 : 35);
			minInsert0=Tools.min(cap, minInsert0);
		}
		
		if(MATE_BY_OVERLAP && !useNormalMode && !useRatioMode){
			System.err.println("\n*** WARNING! Both normal and ratio mode were disabled; using normal mode. ***\n");
			useNormalMode=true;
		}
		
//		assert(false) : useNormalMode+", "+useRatioMode;
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			qtrim=((qtrimLeft||qtrimRight)&&trimq>=0);
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			minReadLength=parser.minReadLength;
			maxReadLength=(parser.maxReadLength<0 ? Integer.MAX_VALUE : parser.maxReadLength);
			untrim=parser.untrim;
		}
		parseCustom=FASTQ.PARSE_CUSTOM;
		if(verbose){
//			assert(false) : "verbose flag is static final; recompile to change it.";
//			BBMergeOverlapper.verbose=true;
		}
		
		if(trimByOverlap){
			join=false;
		}
		
		if(!mm0set){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
		}
		
		if(MAX_MISMATCHES0<MAX_MISMATCHES){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
			System.err.println("MAX_MISMATCHES0 was set to "+MAX_MISMATCHES0+" to remain >=MAX_MISMATCHES");
		}
		
		if(MISMATCH_MARGIN>MAX_MISMATCHES){
			MISMATCH_MARGIN=MAX_MISMATCHES;
			System.err.println("MISMATCH_MARGIN was set to "+MISMATCH_MARGIN+" to remain >=MAX_MISMATCHES");
		}
		
		if(recalibrateQuality){CalcTrueQuality.initializeMatrices();}
		
		if(findAdapterSequence){
			for(int i=0; i<adapterCounts.length; i++){
				for(int j=0; j<adapterCounts[i].length; j++){
					adapterCounts[i][j]=new LongList(150);
				}
			}
		}
		
		if(in2==null && in1!=null && in1.contains("#") && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		
		if(out2==null && out1!=null && out1.contains("#")){
			out2=out1.replaceFirst("#", "2");
			out1=out1.replaceFirst("#", "1");
		}
		
		if(outb2==null && outb1!=null && outb1.contains("#")){
			outb2=outb1.replaceFirst("#", "2");
			outb1=outb1.replaceFirst("#", "1");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outb1, outb2, outinsert, outhist, outhist2, outhist3)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out1+", "+out2+", "+outb1+", "+outb2+", "+outinsert+", "+outhist+", "+outhist2+", "+outhist3+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outb1, outb2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		if(THREADS<1){THREADS=Shared.threads();}
		
		if(useKFilter){
			
			long memory=Runtime.getRuntime().maxMemory();
			long tmemory=Runtime.getRuntime().totalMemory();
			
			final long usable=(long)Tools.max(((memory-100000000-THREADS*8000000)*.70), memory*0.40);
			long mem=usable;
			
//			assert(false) : memory+", "+tmemory+", "+usable;
			
			filterCells=(mem*8)/filterBits;
			if(filterCells<80000000){
				System.err.println("Disabling filter due to low memory.  To enable the filter, increase the -Xmx flag.");
				useKFilter=false;
			}
		}
		useMEEfilter=maxExpectedErrors>0;
		
		Read.VALIDATE_IN_CONSTRUCTOR=(THREADS<16);
	}
	
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
		
		if(useKFilter){
			KmerCountAbstract.CANONICAL=true;
//			assert(false) : filterCells;
			filter=KmerCount7MTA.makeKca(in1, in2, extraFiles, filterK, filterBits, 0, filterCells, filterHashes, filterMinq, true, false, filterReads, 1, 1, 1, 1, null);
			System.err.println("Made prefilter:   \t"+filter.toShortString(filterHashes));
			double uf=filter.usedFraction();
			if(uf>0.7){
				System.err.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.98 ? "incredibly" : uf>0.95 ? "extremely" : uf>0.9 ? "very" : 
					uf>0.8 ? "fairly" : "somewhat")+" full, which may increase false-positive kmers.  Ideal load is under 70% used." +
					"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; " +
						"or increase the values of the minprob flag to reduce spurious kmers.");
			}
		}
		
		runPhase(join, maxReads, false);
		
		double stdev=0;
		if(histTotal!=null){
			stdev=Tools.standardDeviationHistogram(histTotal);
		}
		
		if(outhist!=null){
			StringBuilder sb=new StringBuilder();
			
			sb.append("#Mean\t"+String.format("%.3f", Tools.averageHistogram(histTotal))+"\n");
			sb.append("#Median\t"+Tools.percentile(histTotal, 0.5)+"\n");
			sb.append("#Mode\t"+Tools.calcMode(histTotal)+"\n");
			sb.append("#STDev\t"+String.format("%.3f", Tools.standardDeviationHistogram(histTotal))+"\n");
			sb.append("#InsertSize\tCount\n");
			for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				if(x>0 || !NONZERO_ONLY){
					sb.append(i+"\t"+x+"\n");
				}
			}
			ReadWrite.writeStringInThread(sb, outhist);
		}
		
		if(outhist2!=null){
			StringBuilder sb=new StringBuilder();
			sb.append("#InsertSize\tCount\n");
			int start=(histTotal.length==0 || histTotal[0]>0 ? 0 : 1);
			for(int i=start; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(x+"\n");
			}
			ReadWrite.writeStringInThread(sb, outhist2);
		}

		if(outhist3!=null){
			
			if(!new File(outhist3).exists()){
				StringBuilder sb=new StringBuilder();
				for(int i=0; i<histTotal.length; i+=bin){
					sb.append(i+"\n");
				}
				ReadWrite.writeString(sb, outhist3);
			}
			
			StringBuilder sb=new StringBuilder();
			TextFile tf=new TextFile(outhist3, false, false);
			for(int i=0; i<histTotal.length; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length && i<=insertMaxTotal; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(tf.readLine()+"\t"+x+"\n");
			}
			tf.close();
			ReadWrite.writeStringInThread(sb, outhist3);
		}
		
		if(outAdapter!=null){
			assert(findAdapterSequence);
			writeAdapterConsensus(outAdapter, adapterCounts);
		}
		
		ttotal.stop();
		System.err.println("Total time: "+ttotal+"\n");
		
		long sum=correctCountTotal+incorrectCountTotal;

		double div=100d/readsProcessedTotal;
		double div2=100d/sum;
		System.err.println("Pairs:       \t"+readsProcessedTotal);
		System.err.println("Joined:      \t"+sum+String.format((sum<10000 ? "       " : "   ")+"\t%.3f%%", sum*div));
		System.err.println("Ambiguous:   \t"+ambiguousCountTotal+String.format((ambiguousCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", ambiguousCountTotal*div));
		System.err.println("No Solution: \t"+noSolutionCountTotal+String.format((noSolutionCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", noSolutionCountTotal*div));
		if(minInsert>0){System.err.println("Too Short:   \t"+tooShortCountTotal+String.format((tooShortCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooShortCountTotal*div));}
		if(maxReadLength<Integer.MAX_VALUE){System.err.println("Too Long:    \t"+tooLongCountTotal+String.format((tooLongCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooLongCountTotal*div));}
		
		if(parseCustom){
			System.err.println();
			System.err.println("Correct:     \t"+correctCountTotal+String.format((correctCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", correctCountTotal*div)+String.format("   \t%.3f%% of merged", correctCountTotal*div2));
			System.err.println("Incorrect:   \t"+incorrectCountTotal+String.format((incorrectCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", incorrectCountTotal*div)+String.format("   \t%.3f%% of merged", incorrectCountTotal*div2));
			double snr=Tools.max(correctCountTotal, 0.001)/(Tools.max(incorrectCountTotal, 0.001));
			double snrDB=Tools.mid(-20, 80, 10*Math.log10(snr));
			System.err.println("SNR:         \t"+String.format("%.3f dB", snrDB));
			System.err.println();
			System.err.println("Avg Insert Correct:  \t\t"+String.format("%.1f", (insertSumCorrectTotal)*1d/(correctCountTotal)));
			System.err.println("Avg Insert Incorrect:\t\t"+String.format("%.1f", (insertSumIncorrectTotal)*1d/(incorrectCountTotal)));
		}
		
		System.err.println("\nAvg Insert:          \t\t"+String.format("%.1f", (insertSumCorrectTotal+insertSumIncorrectTotal)*1d/(correctCountTotal+incorrectCountTotal)));
		System.err.println("Standard Deviation:  \t\t"+String.format("%.1f", stdev));
		System.err.println("Mode:                \t\t"+Tools.calcMode(histTotal));
		
		System.err.println();
		System.err.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
		System.err.println("90th percentile:     \t"+Tools.percentile(histTotal, .9));
		System.err.println("75th percentile:     \t"+Tools.percentile(histTotal, .75));
		System.err.println("50th percentile:     \t"+Tools.percentile(histTotal, .5));
		System.err.println("25th percentile:     \t"+Tools.percentile(histTotal, .25));
		System.err.println("10th percentile:     \t"+Tools.percentile(histTotal, .1));
	}
	
	public static void writeAdapterConsensus(String fname, LongList[][] matrix){
		StringBuilder sb=new StringBuilder();
		{
			sb.append(">Read1_adapter\n");
			StringBuilder adapter=new StringBuilder();
			LongList[] lists=matrix[0];
			long max=0;
			int lastBase=-1;
			for(int i=0; true; i++){
				long a=lists[0].get(i);
				long c=lists[1].get(i);
				long g=lists[2].get(i);
				long t=lists[3].get(i);
				long sum=(a+c+g+t);
				max=Tools.max(max, sum);
				if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
				long thresh=(max>100 ? 4+(sum*2)/3 : (sum*2)/3);
				if(a>thresh){
					adapter.append('A');
					lastBase=i;
				}else if(c>thresh){
					adapter.append('C');
					lastBase=i;
				}else if(g>thresh){
					adapter.append('G');
					lastBase=i;
				}else if(t>thresh){
					adapter.append('T');
					lastBase=i;
				}else{
					adapter.append('N');
				}
			}
			if(lastBase<0){sb.append('N');}
			else{
				for(int i=0; i<=lastBase; i++){
					sb.append(adapter.charAt(i));
				}
			}
			sb.append('\n');
		}
		if(matrix.length>1){
			sb.append(">Read2_adapter\n");
			StringBuilder adapter=new StringBuilder();
			LongList[] lists=matrix[1];
			long max=0;
			int lastBase=-1;
			for(int i=0; true; i++){
				long a=lists[0].get(i);
				long c=lists[1].get(i);
				long g=lists[2].get(i);
				long t=lists[3].get(i);
				long sum=(a+c+g+t);
				max=Tools.max(max, sum);
				if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
				long thresh=(max>100 ? 5+(sum*2)/3 : (sum*2)/3);
				if(a>thresh){
					adapter.append('A');
					lastBase=i;
				}else if(c>thresh){
					adapter.append('C');
					lastBase=i;
				}else if(g>thresh){
					adapter.append('G');
					lastBase=i;
				}else if(t>thresh){
					adapter.append('T');
					lastBase=i;
				}else{
					adapter.append('N');
				}
			}
			if(lastBase<0){sb.append('N');}
			else{
				for(int i=0; i<=lastBase; i++){
					sb.append(adapter.charAt(i));
				}
			}
			sb.append('\n');
		}
		ReadWrite.writeString(sb, fname);
	}
	
	public void runPhase(boolean join, long maxReads, boolean perfectonly){
		
		Timer talign=new Timer();
		
		ConcurrentReadOutputStream rosgood=null;
		ConcurrentReadOutputStream rosbad=null;
		ConcurrentReadOutputStream rosinsert=null;
		
		if(out1!=null){
			if(join==true){
				if(out2==null){System.err.println("Writing mergable reads merged.");}
				else{
					System.err.println("WARNING: 2 output files specified even though 'merge=true'.  out2 will be ignored.");
					out2=null;
				}
			}else{
				if(out2==null){System.err.println("Writing mergable reads interleaved.");}
				else{System.err.println("Writing mergable reads unmerged in two files.");}
			}
			
			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosgood=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosgood.start();
		}
		
		if(outb1!=null){

			final FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			final FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosbad=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final int buff=Tools.max(16, 2*THREADS);
			
			String out1=outinsert.replaceFirst("#", "1");

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			
			ReadStreamWriter.HEADER=header();
			final FileFormat ff=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, ordered);
			rosinsert=ConcurrentReadOutputStream.getStream(ff, null, null, null, buff, null, false);
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			System.err.println("Started output threads.");
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		assert(paired);//Fails on empty files.
		if(verbose){System.err.println("Paired: "+paired);}
		
		talign.start();
		
		
		MateThread[] pta=new MateThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new MateThread(cris, rosgood, rosbad, rosinsert, join, perfectonly, trimByOverlap);
			pta[i].start();
		}

		insertMinTotal=999999999;
		insertMaxTotal=0;
		
		readsProcessedTotal=0;
		matedCountTotal=0;
		correctCountTotal=0;
		ambiguousCountTotal=0;
		tooShortCountTotal=0;
		tooLongCountTotal=0;
		incorrectCountTotal=0;
		noSolutionCountTotal=0;
		insertSumCorrectTotal=0;
		insertSumIncorrectTotal=0;
		
		Arrays.fill(histTotal, 0);
		
		for(int i=0; i<pta.length; i++){
			MateThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				readsProcessedTotal+=ct.pairsProcessed;
				matedCountTotal+=ct.matedCount;
				correctCountTotal+=ct.correctCount;
				ambiguousCountTotal+=ct.ambiguousCount;
				tooShortCountTotal+=ct.tooShortCount;
				tooLongCountTotal+=ct.tooLongCount;
				incorrectCountTotal+=ct.incorrectCount;
				noSolutionCountTotal+=ct.noSolutionCount;
				insertSumCorrectTotal+=ct.insertSumCorrect;
				insertSumIncorrectTotal+=ct.insertSumIncorrect;
				
				basesTrimmedTotal+=ct.basesTrimmedT;
				readsTrimmedTotal+=ct.readsTrimmedT;

				insertMinTotal=Tools.min(ct.insertMin, insertMinTotal);
				insertMaxTotal=Tools.max(ct.insertMax, insertMaxTotal);
				
//				System.err.println(ct.insertMin+", "+ct.insertMax);
				
				if(ct.hist!=null){
					for(int h=0; h<ct.hist.length; h++){
						histTotal[h]+=ct.hist[h];
					}
				}
				
				if(findAdapterSequence){
					LongList[][] adapterCountsT=ct.adapterCountsT;
					for(int x=0; x<adapterCounts.length; x++){
						for(int y=0; y<adapterCounts[x].length; y++){
							adapterCounts[x][y].add(adapterCountsT[x][y]);
						}
					}
				}
			}
		}
		
		System.err.println("Finished reading");
		errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad, rosinsert);
		
		talign.stop();
//		System.err.println("Align time: "+talign);
	}
	
	public static final float mergeableFraction(String fname1, String fname2, long numReads, float samplerate){
		long[] hist=makeInsertHistogram(fname1, fname2, numReads, samplerate);
		if(hist==null || hist.length<2){return 0;}
		long sum=Tools.sum(hist);
		return sum<1 ? 0 : (sum-hist[0])/(float)sum;
	}
	
	public static final long[] makeInsertHistogram(String fname1, String fname2, long numReads, float samplerate){
		assert(fname1!=null);
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			if(ff1.stdio()){return null;}
			assert(!ff1.stdio()) : "Standard in is not allowed as input when calculating insert size distributions for files.";
			cris=ConcurrentReadInputStream.getReadInputStream(numReads, true, ff1, ff2);
			cris.setSampleRate(samplerate, 1);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
			if(!cris.paired()){
				ReadWrite.closeStreams(cris);
				return null;
			}
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert(r.mate!=null);
		}

		LongList ll=new LongList(500);
		while(reads!=null && reads.size()>0){

			for(Read r1 : reads){
				int x=findOverlapLoose(r1, r1.mate, false);
				if(x>0){ll.increment(x, 1);}
				else{ll.increment(0, 1);}
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln.id, ln.list.isEmpty());
		ReadWrite.closeStreams(cris);
		return ll.toArray();
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapStrict(final Read r1, final Read r2, boolean ecc){
		final float maxRatio=0.06f;
		final float ratioMargin=10f;
		final float ratioOffset=0.5f;
		
		final float efilterRatio=2f;
		final float efilterOffset=0.45f;
		final float pfilterRatio=0.008f;

		final int minOverlap=8;
		final int minOverlap0=4;
		final int minInsert=50;
		final int minInsert0=35;
		final int entropy=42;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapLoose(final Read r1, final Read r2, boolean ecc){
		
		final float maxRatio=0.12f;
		final float ratioMargin=3f;
		final float ratioOffset=0.45f;
		
		final float efilterRatio=7.5f;
		final float efilterOffset=0.55f;
		final float pfilterRatio=0.000004f;

		final int minOverlap=5;
		final int minOverlap0=6;
		final int minInsert=16;
		final int minInsert0=16;
		final int entropy=28;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}
	
	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlap(final Read r1, final Read r2, final boolean ecc,
			int minOverlap, final int minOverlap0, final int minInsert, final int minInsert0, final int entropy,
			final float maxRatio, final float ratioMargin, final float ratioOffset,
			final float efilterRatio, final float efilterOffset, final float pfilterRatio){
		
		assert(r1!=null && r2!=null);
		if(!r1.validated()){r1.validate(true);}
		if(!r2.validated()){r2.validate(true);}
		
		final boolean swapped;
		if(r2.length()<r1.length()){
			swapped=true;
			r1.swapBasesWithMate();
		}else{
			swapped=false;
		}

		final int len1=r1.length(), len2=r2.length();
		final int minlen=Tools.min(len1, len2);
		
		if(minlen<MIN_OVERLAPPING_BASES || minlen<minInsert){
			if(swapped){r1.swapBasesWithMate();}
			return -1;
		}
		
		int[] rvector=localRvector.get();
		if(rvector==null){
			rvector=new int[5];
			localRvector.set(rvector);
		}
		
		r2.reverseComplement();
		
		int bestInsert=-1;
		int bestBad=999999;
		boolean ambig, tooShort=false;
		
		if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && r1.mapped() && r2.mapped()){
			bestBad=0;
			bestInsert=Read.insertSizeMapped(r1, r2, ignoreMappingStrand);
			ambig=false;
		}else{
			if(entropy>0){
				int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, 3, null, entropy);
				int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, 3, null, entropy);
				minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
			}else{minOverlap=MIN_OVERLAPPING_BASES;}
			if(verbose){System.err.println("minOverlap: "+minOverlap);}
			
			rvector[4]=0;

			int x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, null, null, rvector, minOverlap0, minOverlap, 
					minInsert0, minInsert, maxRatio, ratioMargin, ratioOffset, 0.95f, 0.95f, false);
			bestInsert=x;
			bestBad=rvector[2];
			ambig=(x>-1 ? rvector[4]==1 : false);
		}

		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
		if(bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
			if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
			if(verbose){System.err.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}

		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
		if(pfilterRatio>0 && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
			if(probability<pfilterRatio){bestInsert=-1;}
			if(verbose){System.err.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}

		tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
		
		if(ecc && bestInsert>-1 && !ambig && !tooShort){
			errorCorrectWithInsert(r1, r2, bestInsert);
		}
		
		if(r2!=null){r2.reverseComplement();}
		if(swapped){r1.swapBasesWithMate();}
		if(!ambig && bestInsert>-1){r1.setInsert(bestInsert);}
		
		return ambig ? -1 : bestInsert;
	}
	
	public static int errorCorrectWithInsert(Read r1, Read r2, int insert){
		assert(insert>0);
		int errors=0;
		Read joined=r1.joinRead(insert);
		
		if(joined!=null && joined.length()>0){
			final int lenj=joined.length();
			final int lim1=Tools.min(joined.length(), r1.length());
			final int lim2=lenj-Tools.min(r2.length(), lenj);

			r1.bases=Arrays.copyOfRange(joined.bases, 0, lim1);
			r1.quality=(r1.quality==null ? null : Arrays.copyOfRange(joined.quality, 0, lim1));

			r2.bases=Arrays.copyOfRange(joined.bases, lim2, lenj);
			r2.quality=(r2.quality==null ? null : Arrays.copyOfRange(joined.quality, lim2, lenj));
		}
		return errors;
	}

	public static String header(){
		return "#id\tnumericID\tinsert\tstatus\thashHits\thashMisses\tscore\tsum\tvotes\n";
	}
	
	
	private class MateThread extends Thread{
		
		
		public MateThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosgood_, ConcurrentReadOutputStream rosbad_, ConcurrentReadOutputStream rosinsert_,
				boolean joinReads_, boolean joinperfectonly_, boolean trimByOverlap_) {
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosinsert=rosinsert_;
			joinReads=joinReads_;
			trimReadsByOverlap=trimByOverlap_;
			joinperfectonly=joinperfectonly_;
			
			if(useEntropy){
				kmerCounts=new short[1<<(2*entropyK)];
			}else{
				kmerCounts=null;
			}
			
			if(findAdapterSequence){
				for(int i=0; i<adapterCountsT.length; i++){
					for(int j=0; j<adapterCountsT[i].length; j++){
						adapterCountsT[i][j]=new LongList(150);
					}
				}
			}
		}
		
		
		@Override
		public void run(){
			processReads();
		}

		private void processReads() {
			assert(USE_MAPPING || MATE_BY_OVERLAP);

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null);
			}


			while(reads!=null && reads.size()>0){

				ArrayList<Read> listg=(rosgood==null ? null : new ArrayList<Read>(reads.size()));
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>(reads.size()));
				ArrayList<Read> listi=(rosinsert==null ? null : new ArrayList<Read>(reads.size()));

				for(Read r1 : reads){
					processReadPair(r1, r1.mate, listg, listb, listi);
				}

				if(rosgood!=null){rosgood.add(listg, ln.id);}
				if(rosbad!=null){rosbad.add(listb, ln.id);}
				if(rosinsert!=null){rosinsert.add(listi, ln.id);}

				//			System.err.println("returning list");
				cris.returnList(ln.id, ln.list.isEmpty());
				//			System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			System.err.println("reads: "+(reads==null ? "null" : reads.size()));
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}
		
		
		private final int processReadPair(final Read r1, final Read r2, ArrayList<Read> listg, ArrayList<Read> listb, ArrayList<Read> listi){
			
			assert(r1!=null);
			if(!r1.validated()){r1.validate(true);}
			if(r2==null){
				listb.add(r1);
				return r1.length();
			}
			if(!r2.validated()){r2.validate(true);}
			
			if(iupacToN){
				if(r1!=null){r1.convertUndefinedTo((byte)'N');}
				if(r2!=null){r2.convertUndefinedTo((byte)'N');}
			}
			
			if(recalibrateQuality){
				CalcTrueQuality.recalibrate(r1);
				CalcTrueQuality.recalibrate(r2);
			}
			
			final boolean swapped;
			if(useRatioMode && r2.length()<r1.length()){
				swapped=true;
				r1.swapBasesWithMate();
			}else{
				swapped=false;
			}
			
			pairsProcessed++;
			
			TrimRead tr1=null, tr2=null;
			if(qtrim){
				if(untrim){
					tr1=TrimRead.trim(r1, qtrimLeft, qtrimRight, trimq, 1);
					int x1=(tr1==null ? 0 : tr1.leftTrimmed+tr1.rightTrimmed);
					basesTrimmedT+=x1;
					readsTrimmedT+=(x1>0 ? 1 : 0);

					tr2=TrimRead.trim(r2, qtrimLeft, qtrimRight, trimq, 1);
					int x2=(tr2==null ? 0 : tr2.leftTrimmed+tr2.rightTrimmed);
					basesTrimmedT+=x2;
					readsTrimmedT+=(x2>0 ? 1 : 0);
				}else{
					int x1=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
					basesTrimmedT+=x1;
					readsTrimmedT+=(x1>0 ? 1 : 0);

					int x2=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
					basesTrimmedT+=x2;
					readsTrimmedT+=(x2>0 ? 1 : 0);
				}
			}
			

			final int len1=r1.length(), len2=r2.length();
			
			if(minReadLength>0 && len1<minReadLength && len2<minReadLength){
				basesTrimmedT+=(len1+len2);
				if(swapped){r1.swapBasesWithMate();}
				return -1;
			}
			
			if(r1.quality!=null || r2.quality!=null){
				if(minAvgQuality>0){
					if(r1.avgQuality(false, minAvgQualityBases)<minAvgQuality || r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){
						basesTrimmedT+=(len1+len2);
						if(swapped){r1.swapBasesWithMate();}
						return -1;
					}
				}
				if(useMEEfilter && useQuality){
					int maxBasesToConsider=Tools.min(Tools.max(len1, len2), len1+len2-minInsert);
					if(r1.expectedTipErrors(false, maxBasesToConsider)>maxExpectedErrors || r2.expectedTipErrors(false, maxBasesToConsider)>maxExpectedErrors){
						basesTrimmedT+=(len1+len2);
						if(swapped){r1.swapBasesWithMate();}
						return -1;
					}
				}
			}
			
			int trueSize=-1;
			
			if(parseCustom){
				if(r1.id.startsWith("insert=")){
					trueSize=GradeMergedReads.parseInsert(r1.id);
				}else{
					r1.setMapped(true);
					r2.setMapped(true);
					trueSize=Read.insertSizeMapped(r1, r2, ignoreMappingStrand);
				}
				if(verbose){System.err.println("True Insert: "+trueSize);}
			}
			
			r2.reverseComplement();
			
			int bestInsert=-1;
			int bestBad=999999;
			boolean ambig, tooShort=false, tooLong=false;
			
			final byte[] qual1=r1.quality, qual2=r2.quality;
			if(!useQuality){//strip qualities
				r1.quality=r2.quality=null;
			}
			
			if(len2<1){
				bestBad=0;
				bestInsert=len1;
				assert(len1==r1.insert()) : len1+" != "+r1.insert()+"; actual = "+trueSize;
				ambig=false;
			}else if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && ((r1.mapped() || r1.synthetic()) && (r2.mapped() || r2.synthetic()))){
				bestBad=0;
				bestInsert=trueSize;
				ambig=false;
			}else if(SKIP_MATED_READS && r1.insertvalid() && r1.insert()>0){
				bestBad=0;
				bestInsert=r1.insert();
				ambig=false;
			}else{
				if(MATE_BY_OVERLAP){
					
					boolean ambigRM;
					int bestBadRM, bestInsertRM;
					
					final int minOverlap;
					if(useEntropy){
						if(loose){
							int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, entropyK, kmerCounts, minEntropyScore);
							int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, entropyK, kmerCounts, minEntropyScore);
							float errorRate=r1.expectedErrors(false, len1)/len1+r2.expectedErrors(false, len2)/len2;
							minOverlap=(int)(Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b))*(1+Tools.min(0.04f, errorRate)*4f));
						}else{
							int a=BBMergeOverlapper.calcMinOverlapByEntropyTail(r1.bases, entropyK, kmerCounts, minEntropyScore);
							int b=BBMergeOverlapper.calcMinOverlapByEntropyHead(r2.bases, entropyK, kmerCounts, minEntropyScore);
							minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
						}
					}else{minOverlap=MIN_OVERLAPPING_BASES;}
					if(verbose){System.err.println("minOverlap: "+minOverlap);}
					
					boolean ambigNM=false;
					int bestInsertNM=-1;
					int bestBadNM=999999;
					
					if(lowercaseAdapters && Character.isLowerCase(r1.bases[r1.length()-1]) && Character.isLowerCase(r2.bases[0])){
						final int lower1=r1.trailingLowerCase(), lower2=r2.leadingLowerCase();

						final int upper1=r1.length()-lower1, upper2=r2.length()-lower2;
						final int newlen=Tools.min(upper1, upper2);
						int good=0, bad=0;

						for(int i=0; i<newlen; i++){
							byte a=r1.bases[i];
							byte b=r2.bases[i+lower2];
							if(a!='N' && b!='N'){
								if(a==b){good++;}
								else{bad++;}
							}
						}
						if(bad*4<=good){
							bestInsertNM=newlen;
							ambigNM=false;
							bestBadNM=bad;
						}
					}
					
					if(aprob==null || aprob.length<len1){aprob=new float[len1];}
					if(bprob==null || bprob.length<len2){bprob=new float[len2];}
					
					if(useRatioMode && bestInsertNM<0){
						int min0=MIN_OVERLAPPING_BASES_0-MIN_OVERLAPPING_BASES_RATIO_REDUCTION;
						int min=minOverlap-MIN_OVERLAPPING_BASES_RATIO_REDUCTION;
						int x=-1;
						rvector[4]=0;
						
//						float errorRate=Tools.max(r1.expectedErrors(false, len1)/len1, r2.expectedErrors(false, len2)/len2); 
//						float ratioMargin=(false ? RATIO_MARGIN*(1+Tools.min(0.15f, errorRate*1.2f)) : RATIO_MARGIN);
//						float maxRatio=(false ? MAX_RATIO*(1+Tools.min(0.15f, errorRate*2f)) : MAX_RATIO);
						
						float ratioMargin=RATIO_MARGIN;
						float maxRatio=MAX_RATIO;
						
						boolean overlapped=false;
						if(overlapUsingQuality && r1.quality!=null && r2.quality!=null){
							overlapped=true;
							x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, 
									min0, min, minInsert0, minInsert, maxRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, true);
							if(verbose){System.err.println("Result from ratiomode1:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
						}
						if(!overlapped || (overlapWithoutQuality && (x<0 || rvector[4]==1))){
							x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, min0, min, 
									minInsert0, minInsert, maxRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, false);
							if(verbose){System.err.println("Result from ratiomode2:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
						}
						bestInsertRM=x;
						bestBadRM=rvector[2];
						ambigRM=(x>-1 ? rvector[4]==1 : false);
						
//						if(TRIM_ON_OVERLAP_FAILURE>0 && (bestInsertRM<0 || ambigRM)){
//							Serializable old1=r1.obj;
//							Serializable old2=r2.obj;
//							tr1=TrimRead.trim(r1, false, true, 10, 1+len1*4/10); //r1.length());
//							tr2=TrimRead.trim(r2, true, false, 10, 1+len2*4/10); //r2.length());
//							r1.obj=old1;
//							r2.obj=old2;
//							overlapped=false;
//							if(overlapUsingQuality && r1.quality!=null && r2.quality!=null){
//								overlapped=true;
//								x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, 
//										min0, min, minInsert0, minInsert, maxRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, true);
//								if(verbose){System.err.println("Result from ratiomode1:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
//							}
//							if(!overlapped || (overlapWithoutQuality && (x<0 || rvector[4]==1))){
//								x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, min0, min, 
//										minInsert0, minInsert, maxRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, false);
//								if(verbose){System.err.println("Result from ratiomode2:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
//							}
//							
//							if(x>-1 && rvector[4]!=1){
//								bestInsertRM=x;
//								bestBadRM=rvector[2];
//								ambigRM=(x>-1 ? rvector[4]==1 : false);
//							}else{
//								if(tr1!=null){tr1.untrim();}
//								if(tr2!=null){tr2.untrim();}
//							}
//						}
						
					}else{
						bestInsertRM=-1;
						bestBadRM=0;
						ambigRM=false;
					}
					
					if(useNormalMode && ((!requireRatioMatch && (bestInsertRM<0 || ambigRM)) || (requireRatioMatch && (bestInsertRM>0 && !ambigRM)))){
						bestInsertNM=-1;
						assert(QUAL_ITERS>0);
						final int maxQualIters=(r1.quality==null || r2.quality==null ? 1 : QUAL_ITERS);
						final int maxTrims=(r1.quality==null || r2.quality==null ? 0 : TRIM_ON_OVERLAP_FAILURE);

						for(int i=0; i<maxQualIters && bestInsertNM<0 /*&& !ambigNM*/; i++){

							int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-i, minOverlap+i, 
									minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, (byte)(MIN_QUALITY-2*i));
							if(x>-1){
								bestInsertNM=x;
								bestBadNM=rvector[2];
								ambigNM=(rvector[4]==1);
								break;
							}
						}


						if(loose && bestInsertNM<0){//TODO check for estimated number of overlap errors
							int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap+2, 
									minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0+1, MAX_MISMATCHES+1, MIN_QUALITY-1);
							if(x>-1){
								bestInsertNM=x;
								bestBadNM=rvector[2];
								ambigNM=(rvector[4]==1);
							}
						}

						for(int trims=0, q=trimq; trims<maxTrims && !qtrim && bestInsertNM<0 /*&& !ambigNM*/; trims++, q+=8){
							Serializable old1=r1.obj;
							Serializable old2=r2.obj;
							tr1=TrimRead.trim(r1, false, true, q, 1+len1*4/10); //r1.length());
							tr2=TrimRead.trim(r2, true, false, q, 1+len2*4/10); //r2.length());
							r1.obj=old1;
							r2.obj=old2;
							if(tr1!=null || tr2!=null){
								int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap, 
										minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, MIN_QUALITY);
								if(x>-1){
									bestInsertNM=x;
									bestBadNM=rvector[2];
									ambigNM=(rvector[4]==1);
									trims=maxTrims;
								}else{
									if(tr1!=null){tr1.untrim();}
									if(tr2!=null){tr2.untrim();}
								}
							}
						}
						if(verbose){System.err.println("Result from normalmode:  \tinsert="+bestInsertNM+", bad="+bestBadNM+", ambig="+ambigNM);}
					}
					
					if(requireRatioMatch && useNormalMode && useRatioMode){
						ambig=ambigRM || ambigNM;
						bestBad=bestBadRM;
						bestInsert=(bestInsertNM==bestInsertRM ? bestInsertNM : -1);

						if(verbose){System.err.println("Result after rrm:  \tinsert="+bestInsertNM+", bad="+bestBadNM+", ambig="+ambigNM);}
					}else if(useRatioMode && bestInsertRM>-1 && !ambigRM){
						ambig=ambigRM;
						bestBad=bestBadRM;
						bestInsert=bestInsertRM;
					}else{
						ambig=ambigNM;
						bestBad=bestBadNM;
						bestInsert=bestInsertNM;
					}

					//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
					if(useEfilter && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
						float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
						if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
						if(verbose){System.err.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
					}

					//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
					if(pfilterRatio>0 && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
						float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
						if(probability<pfilterRatio){bestInsert=-1;}
						if(verbose){System.err.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
					}
					
//					assert(ambig || bestInsert<0 || trueSize==bestInsert) : bestInsert+", "+trueSize+", "+r1.numericID+"\n"+r1.toFastq()+"\n"+r2.toFastq()+"\n";
					
				}else{
					ambig=false;
					bestInsert=-1;
				}
			}
			
			if(!useQuality){//restore qualities
				r1.quality=qual1;
				r2.quality=qual2;
			}

			tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
			tooLong=(!ambig && !tooShort && bestInsert>0 && bestInsert>maxReadLength);

			if(joinperfectonly && bestBad>0){ambig=true;}

			Read joined=null;
			
			if(bestInsert>-1 && !ambig && !tooShort){
				if(joinReads || useKFilter){
					joined=r1.joinRead(bestInsert);
					if(useKFilter){
						int cov=BBMergeOverlapper.minCoverage(joined, filter, filterK, true, filterCutoff);
						if(cov<filterCutoff){bestInsert=-1;}
						if(verbose){System.err.println("Result after kfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
					}
					if(bestInsert<0){joined=null;}
					if(ecc && joined!=null && joined.length()>0){
						final int lenj=joined.length();
						final int lim1=Tools.min(joined.length(), r1.length());
						final int lim2=lenj-Tools.min(len2, lenj);

						//					assert(len1==Arrays.copyOfRange(joined.bases, 0, lim1).length || Arrays.copyOfRange(joined.bases, 0, lim1).length==bestInsert && bestInsert<len1) : r1.length()+", "+lim1;

						r1.bases=Arrays.copyOfRange(joined.bases, 0, lim1);
						r1.quality=(r1.quality==null ? null : Arrays.copyOfRange(joined.quality, 0, lim1));


						//					assert(r2.length()==Arrays.copyOfRange(joined.bases, lim2, lenj).length || Arrays.copyOfRange(joined.bases, lim2, lenj).length==bestInsert && bestInsert<len2);

						r2.bases=Arrays.copyOfRange(joined.bases, lim2, lenj);
						r2.quality=(r2.quality==null ? null : Arrays.copyOfRange(joined.quality, lim2, lenj));

						//					if(r1.length()<150){
						//						System.err.println("\n\n"+r1.length()+"\n"+r1.toFastq()+"\n"+r2.toFastq()+"\n");
						//					}

						//					assert(r1.length()==r2.length()) : r1.length()+"\n"+r1.toFastq()+"\n"+r2.toFastq()+"\n";
						//					assert(r1.length()==150) : r1.length()+"\n"+r1.toFastq()+"\n"+r2.toFastq()+"\n";
					}
					if(!joinReads){joined=null;}
				}else if(ecc){
					errorCorrectWithInsert(r1, r2, bestInsert);
				}
			}
			
			if(ambig){ambiguousCount++;}
			else if(tooShort){tooShortCount++;}
			else if(tooLong){tooLongCount++;}
			else if(bestInsert==-1){noSolutionCount++;}
			else if(bestInsert==trueSize){correctCount++;insertSumCorrect+=bestInsert;}
			else{incorrectCount++;insertSumIncorrect+=bestInsert;}

			if(bestInsert>-1 && !ambig){
				insertMin=Tools.min(bestInsert, insertMin);
				insertMax=Tools.max(bestInsert, insertMax);
				hist[Tools.min(bestInsert, hist.length-1)]++;
			}
			r1.setInsert(ambig ? -1 : bestInsert);

			if(OUTPUT_FAILED || bestInsert>-1){

				if(untrim && (ambig || bestInsert<0 || !joinReads)){
					if(tr1!=null){tr1.untrim();}
					if(tr2!=null){tr2.untrim();}
				}

				if((ambig || bestInsert<0 || tooShort) && (listb!=null || !MIX_BAD_AND_GOOD)){
					if(listb!=null){
						if(listb!=null){listb.add(r1);}
					}
				}else{
					if(listg!=null){
						if(joined==null){
							joined=r1;
							if(joinReads){joined=r1.joinRead(bestInsert);}
							else if(trimReadsByOverlap){
								int trimLim=bestInsert-1;
								if(trimLim<r1.length()){
									if(verbose){System.err.println("Overlap right trimming r1 to "+0+", "+(trimLim));}
									int x=TrimRead.trimToPosition(r1, 0, trimLim, 1);
									if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r1.bases));}
								}
								if(trimLim<r2.length()){
									if(verbose){System.err.println("Overlap right trimming r2 to "+0+", "+(trimLim));}
									int x=TrimRead.trimToPosition(r2, 0, trimLim, 1);
									if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r2.bases));}
								}
							}
						}
						listg.add(joined);
					}
				}

				if(listi!=null){
					StringBuilder sb=new StringBuilder(40);
					sb.append(r1.id==null ? r1.numericID+"" : r1.id).append('\t');
					sb.append(r1.numericID).append('\t');

					sb.append(bestInsert);
					sb.append('\t');

					if(bestInsert<0){sb.append('F');}//Failed
					else if(ambig){sb.append('A');} //Ambiguous
					else if(tooShort){sb.append('S');} //Short
					else if(bestInsert>0 && bestBad<1){sb.append('P');} //Perfect
					else{sb.append('I');}//Imperfect

					if(bestInsert>0){
						sb.append("\t"+bestBad);
					}
					r1.obj=sb;
					listi.add(r1);
				}
			}
			if(r2!=null){r2.reverseComplement();}
			
			if(swapped){
				r1.swapBasesWithMate();
				if(joined!=null){joined.reverseComplement();}
			}
			
			if(findAdapterSequence && bestInsert>0 && !ambig){
				storeAdapterSequence(r1, bestInsert);
				storeAdapterSequence(r2, bestInsert);
			}
			
			return bestInsert;
		}
		
		private void storeAdapterSequence(Read r, int insert){
			LongList[] lists=adapterCountsT[r.pairnum()];
			byte[] bases=r.bases;
			for(int i=insert, j=0; i<bases.length; i++, j++){
				byte b=bases[i];
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					lists[num].increment(j);
				}
			}
		}
		
		
		/*--------------------------------------------------------------*/
		
		final LongList[][] adapterCountsT=new LongList[2][4];
		
		private final int[] rvector=new int[5];
		
		final long[] hist=new long[histlen];
		final short[] kmerCounts;
		
		private float[] aprob, bprob;

		long pairsProcessed=0;
		long matedCount=0;
		long correctCount=0;
		long ambiguousCount=0;
		long tooShortCount=0;
		long tooLongCount=0;
		long incorrectCount=0;
		long noSolutionCount=0;
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		int insertMax=0;
		int insertMin=999999999;
		
		long basesTrimmedT=0;
		long readsTrimmedT=0;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream rosgood;
		private final ConcurrentReadOutputStream rosbad;
		private final ConcurrentReadOutputStream rosinsert;
		
		private final boolean joinReads;
		private final boolean joinperfectonly;
		private final boolean trimReadsByOverlap;
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1;
	private String in2;

	private String out1=null;
	private String out2=null;
	private String outb1=null;
	private String outb2=null;
	private String outinsert=null;
	private String outhist=null;
	private String outhist2=null;
	private String outhist3=null;
	private String outAdapter=null;
	
	private long maxReads=-1;
	private boolean join=true;
	private boolean ecc=false;
	private boolean trimByOverlap=false;
	
	private float pfilterRatio=0.00002f;
	private float efilterRatio=6f;
	private float efilterOffset=0.05f;
	private boolean useEfilter=true;
	private boolean useMEEfilter=false;
	
	private boolean ordered=false;
	private boolean overlapUsingQuality=false;
	private boolean overlapWithoutQuality=true;
	private boolean useKFilter=false;
	private int filterBits=1;
	private int filterCutoff=(1<<filterBits)-1;
	private long filterCells=-1;
	private int filterK=31;
	private int filterHashes=1;
	private long filterReads=-1;
	private KCountArray filter;
	private int filterMinq=2;
	private ArrayList<String> extraFiles;

	private boolean useEntropy=true;
	private int entropyK=3;
	private int minEntropyScore=39;//30 loose;//39 normal;//44 strict;
	
	private long sampleseed=-1;
	private float samplerate=1;
	
	private boolean findAdapterSequence=false;
	
	private final LongList[][] adapterCounts=new LongList[2][4];
	
	/*--------------------------------------------------------------*/
	
	private static ThreadLocal<int[]> localRvector=new ThreadLocal<int[]>();
	
	static boolean errorState=false;

	/** Recalibrate quality scores using matrices */
	static boolean recalibrateQuality=false;
	static boolean useQuality=true;
	static boolean qtrimRight=false;
	static boolean qtrimLeft=false;
	static boolean untrim=false;
	static byte trimq=6;
	static byte minAvgQuality=0;
	static int minAvgQualityBases=0;
	static float maxExpectedErrors=0;
	static int minReadLength=0;
	static int maxReadLength=-1;
	static int minInsert=35;
	static int minInsert0=-1;
	static boolean qtrim=false;
	static int TRIM_ON_OVERLAP_FAILURE=1;
	static int QUAL_ITERS=3;
	static boolean parseCustom=false;
	
	static boolean strict=false;
	static boolean vstrict=false;
	static boolean ustrict=false;
	static boolean xstrict=false;
	static boolean loose=false;
	static boolean vloose=false;
	static boolean uloose=false;
	static boolean xloose=false;
	static boolean fast=false;
	
	/** If true, interpret lowercase bases as adapter sequence */
	static boolean lowercaseAdapters=false;
	
	private static final int histlen=2000; 
	static long[] histTotal=new long[histlen];
	static int bin=1;

	static long readsProcessedTotal=0;
	static long matedCountTotal=0;
	static long correctCountTotal=0;
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long tooLongCountTotal=0;
	static long incorrectCountTotal=0;
	static long noSolutionCountTotal=0;
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	static long basesTrimmedTotal=0;
	static long readsTrimmedTotal=0;
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;
	
	private static int MIN_OVERLAPPING_BASES=11;
	private static int MIN_OVERLAPPING_BASES_0=8;
	private static int MISMATCH_MARGIN=2;
	private static int MIN_OVERLAPPING_BASES_RATIO_REDUCTION=3;
	
	static boolean useRatioMode=true;
	static boolean useNormalMode=false;
	static boolean requireRatioMatch=false;
	static float MAX_RATIO=0.09f;
	static float RATIO_MARGIN=5.5f;
	static float RATIO_OFFSET=0.55f;
	
	public static int MAX_MISMATCHES=3;
	public static int MAX_MISMATCHES0=3;
	public static byte MIN_QUALITY=10;
	
	/** Skip alignment and calculate insert from mapping info */ 
	protected static boolean USE_MAPPING=false;
	protected static final boolean ignoreMappingStrand=false;
	
	private static boolean MATE_BY_OVERLAP=true;
	private static boolean SKIP_MATED_READS=false;
	private static boolean OUTPUT_FAILED=true;
	private static boolean MIX_BAD_AND_GOOD=false;
	private static boolean NONZERO_ONLY=true;
	
	private static boolean overwrite=true;
	private static boolean append=false;
	private static final boolean verbose=false;
	
	private static boolean iupacToN=false;
	
	private static int THREADS=-1;
	private static float version=7.3f;
	
}
