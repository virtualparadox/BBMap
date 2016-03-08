package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;
import stream.ReadStreamWriter;

import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;

import align2.ListNum;
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
		BBMerge mr=new BBMerge(args);
		mr.process();
	}
	
	public BBMerge(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		System.err.println("BBMerge version "+version);
		
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
		
		FASTQ.PARSE_CUSTOM=false;
		ReadWrite.MAX_ZIP_THREADS=Shared.THREADS-1;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		
		boolean trimRight_=false;
		boolean trimLeft_=false;
		byte trimq_=trimq;
		int minReadLength_=0;
		
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
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("minoverlappingbases") || a.equals("minoverlapbases") || a.equals("minoverlap")){
				MIN_OVERLAPPING_BASES=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases0") || a.equals("minoverlapbases0") || a.equals("minoverlap0")){
				MIN_OVERLAPPING_BASES_0=Integer.parseInt(b);
			}else if(a.equals("minoverlapinsert") || a.equals("minoi")){
				BBMergeOverlapper.MIN_OVERLAP_INSERT=Integer.parseInt(b);
			}else if(a.equals("minqo") || a.equals("minq")){
				BBMergeOverlapper.MIN_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("maxbadbases") || a.equals("maxbad")){
				BBMergeOverlapper.DEFAULT_BADLIMIT=Integer.parseInt(b);
			}else if(a.equals("usemapping")){
				BBMergeOverlapper.USE_MAPPING=Tools.parseBoolean(b);
			}else if(a.equals("bin")){
				bin=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Shared.THREADS=Integer.parseInt(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
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
			}else if(a.startsWith("outhist") || a.startsWith("hist") || a.equals("ihist")){
				outhist=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist2") || a.equals("hist2")){
				outhist2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist3") || a.equals("hist3")){
				outhist3=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outputfailed")){
				OUTPUT_FAILED=Tools.parseBoolean(b);
			}else if(a.equals("mix")){
				MIX_BAD_AND_GOOD=Tools.parseBoolean(b);
			}else if(a.equals("testinterleaved")){
				FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				System.err.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
			}else if(a.equals("forceinterleaved")){
				FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
				System.err.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("join") || a.equals("merge")){
				join=Tools.parseBoolean(b);
			}else if(a.equals("useoverlap") || a.equals("usebases") || a.equals("matebyoverlap") || a.equals("matebybases")){
				MATE_BY_OVERLAP=Tools.parseBoolean(b);
			}else if(a.startsWith("skipmated")){
				SKIP_MATED_READS=Tools.parseBoolean(b);
			}else if(a.equals("parsecustom")){
				FASTQ.PARSE_CUSTOM=Tools.parseBoolean(b);
				System.err.println("Setting FASTQ.PARSE_CUSTOM to "+FASTQ.PARSE_CUSTOM);
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
			}else if(a.equals("untrim")){
				untrim=Tools.parseBoolean(b);
			}else if(a.equals("trim") || a.equals("qtrim")){
				if(b==null){trimRight_=trimLeft_=true;}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){trimLeft_=true;trimRight_=false;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){trimLeft_=false;trimRight_=true;}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){trimLeft_=trimRight_=true;}
				else{trimRight_=trimLeft_=Tools.parseBoolean(b);}
			}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
				if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
					TrimRead.optimalMode=true;
					TrimRead.optimalBias=Float.parseFloat(b);
					assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
				}else{
					TrimRead.optimalMode=Tools.parseBoolean(b);
				}
			}else if(a.equals("trimright") || a.equals("qtrimright")){
				trimRight_=Tools.parseBoolean(b);
			}else if(a.equals("trimleft") || a.equals("qtrimleft")){
				trimLeft_=Tools.parseBoolean(b);
			}else if(a.equals("trimq") || a.equals("trimquality")){
				trimq_=Byte.parseByte(b);
			}else if(a.equals("q102matrix") || a.equals("q102m")){
				CalcTrueQuality.q102matrix=b;
			}else if(a.equals("qbpmatrix") || a.equals("bqpm")){
				CalcTrueQuality.qbpmatrix=b;
			}else if(a.equals("adjustquality") || a.equals("adjq")){
				TrimRead.ADJUST_QUALITY=Tools.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minReadLength_=Integer.parseInt(b);
			}else if(a.equals("mi") || a.equals("minins") || a.equals("mininsert")){
				minInsert=Integer.parseInt(b);
			}else if(a.equals("ignorebadquality") || a.equals("ibq")){
				FASTQ.IGNORE_BAD_QUALITY=Tools.parseBoolean(b);
			}else if(a.equals("ascii") || a.equals("quality") || a.equals("qual")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;}
				else{x=(byte)Integer.parseInt(b);}
				qin=qout=x;
			}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=true;}
				else{x=(byte)Integer.parseInt(b);}
				qin=x;
			}else if(a.equals("asciiout") || a.equals("qualityout") || a.equals("qualout") || a.equals("qout")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY_OUT=true;}
				else{x=(byte)Integer.parseInt(b);}
				qout=x;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
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
		
		trimRight=trimRight_;
		trimLeft=trimLeft_;
		trimq=trimq_;
		minReadLength=minReadLength_;
		qtrim=trimLeft_||trimRight_;
		
		if(qin!=-1 && qout!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY=false;
		}else if(qin!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.DETECT_QUALITY=false;
		}else if(qout!=-1){
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY_OUT=false;
		}
		
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		if(THREADS<1){THREADS=Shared.THREADS;}
	}
	
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
		
		runPhase(join, maxReads, false);
		
		if(outhist!=null){
			StringBuilder sb=new StringBuilder();
			for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
				int x=0;
				int y=0;
				for(int j=i; j<i+bin && j<histTotal.length; j++){
					x+=histTotal[j];
					y++;
				}
				x=(x+bin-1)/y;
				sb.append(i+"\t"+x+"\n");
			}
			ReadWrite.writeStringInThread(sb, outhist);
		}

		if(outhist2!=null){
			StringBuilder sb=new StringBuilder();
			for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
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
		
		
		ttotal.stop();
		System.err.println("Total time: "+ttotal+"\n");
		
		long sum=correctCountTotal+incorrectCountTotal;
		
		double div=100d/readsProcessedTotal;
		System.err.println("Reads:       \t"+readsProcessedTotal);
		System.err.println("Joined:      \t"+sum+String.format((sum<10000 ? "       " : "   ")+"\t%.3f%%", sum*div));
		if(FASTQ.PARSE_CUSTOM){
			System.err.println("Correct:     \t"+correctCountTotal+String.format((correctCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", correctCountTotal*div));
			System.err.println("Incorrect:   \t"+incorrectCountTotal+String.format((incorrectCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", incorrectCountTotal*div));
		}
		System.err.println("Ambiguous:   \t"+ambiguousCountTotal+String.format((ambiguousCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", ambiguousCountTotal*div));
		System.err.println("No Solution: \t"+noSolutionCountTotal+String.format((noSolutionCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", noSolutionCountTotal*div));
		if(minInsert>0){System.err.println("Too Short:   \t"+tooShortCountTotal+String.format((tooShortCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooShortCountTotal*div));}
		System.err.println("Avg Insert:          \t\t"+String.format("%.1f", (insertSumCorrectTotal+insertSumIncorrectTotal)*1d/(correctCountTotal+incorrectCountTotal)));
		if(FASTQ.PARSE_CUSTOM){
			System.err.println("Avg Insert Correct:  \t\t"+String.format("%.1f", (insertSumCorrectTotal)*1d/(correctCountTotal)));
			System.err.println("Avg Insert Incorrect:\t\t"+String.format("%.1f", (insertSumIncorrectTotal)*1d/(incorrectCountTotal)));
		}
		
		System.err.println();
		System.err.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
		System.err.println("90th percentile:     \t"+Tools.percentile(histTotal, .9));
		System.err.println("50th percentile:     \t"+Tools.percentile(histTotal, .5));
		System.err.println("10th percentile:     \t"+Tools.percentile(histTotal, .1));
	}
	
	public void runPhase(boolean join, long maxReads, boolean perfectonly){
		
		Timer talign=new Timer();
		
		RTextOutputStream3 rosgood=null;
		RTextOutputStream3 rosbad=null;
		RTextOutputStream3 rosinsert=null;
		
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

			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(8, 2*THREADS);
			rosgood=new RTextOutputStream3(ff1, ff2, null, null, buff, null, false);
			rosgood.start();
		}
		
		if(outb1!=null){

			final FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, append, true);
			final FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, append, true);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosbad=new RTextOutputStream3(ff1, ff2, null, null, buff, null, false);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final int buff=Tools.max(16, 2*THREADS);
			
			String out1=outinsert.replaceFirst("#", "1");

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			
			ReadStreamWriter.HEADER=header();
			final FileFormat ff=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, true);
			rosinsert=new RTextOutputStream3(ff, null, null, null, buff, null, false);
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			System.err.println("Started output threads.");
		}
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		assert(paired);
		if(verbose){System.err.println("Paired: "+paired);}
		
		talign.start();
		
		
		MateThread[] pta=new MateThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new MateThread(cris, rosgood, rosbad, rosinsert, join, perfectonly);
			pta[i].start();
		}

		insertMinTotal=999999999;
		insertMaxTotal=0;
		
		readsProcessedTotal=0;
		matedCountTotal=0;
		correctCountTotal=0;
		ambiguousCountTotal=0;
		tooShortCountTotal=0;
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
				
				readsProcessedTotal+=ct.readsProcessed;
				matedCountTotal+=ct.matedCount;
				correctCountTotal+=ct.correctCount;
				ambiguousCountTotal+=ct.ambiguousCount;
				tooShortCountTotal+=ct.tooShortCount;
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
			}
		}
		
		System.err.println("Finished reading");
		errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad, rosinsert);
		
		talign.stop();
		System.err.println("Align time: "+talign);
	}
	


	public static String header(){
		return "#id\tnumericID\tinsert\tstatus\thashHits\thashMisses\tscore\tsum\tvotes\n";
	}
	
	
	private static class MateThread extends Thread{
		
		
		public MateThread(ConcurrentReadStreamInterface cris_, RTextOutputStream3 rosgood_, RTextOutputStream3 rosbad_, RTextOutputStream3 rosinsert_,
				boolean joinReads_, boolean joinperfectonly_) {
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosinsert=rosinsert_;
			joinReads=joinReads_;
			joinperfectonly=joinperfectonly_;
		}
		
		
		@Override
		public void run(){
			processMate();
		}

		private void processMate() {
			
			final boolean USE_MAPPING=BBMergeOverlapper.USE_MAPPING;
			final boolean ignoreMappingStrand=BBMergeOverlapper.IGNORE_MAPPING_STRAND;
			assert(USE_MAPPING || MATE_BY_OVERLAP);

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null);
			}


			while(reads!=null && reads.size()>0){

				ArrayList<Read> listg=(rosgood==null ? null : new ArrayList<Read>());
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>());
				ArrayList<Read> listi=(rosinsert==null ? null : new ArrayList<Read>());

				for(Read r1 : reads){
					final Read r2=r1.mate;
					
					TrimRead tr1=null, tr2=null;
					
					boolean remove=false;
					if(qtrim){
						if(untrim){
							if(r1!=null){
								tr1=TrimRead.trim(r1, trimLeft, trimRight, trimq, 1);
								int x=(tr1==null ? 0 : tr1.leftTrimmed+tr1.rightTrimmed);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2!=null){
								tr2=TrimRead.trim(r2, trimLeft, trimRight, trimq, 1);
								int x=(tr2==null ? 0 : tr2.leftTrimmed+tr2.rightTrimmed);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
						}else{
							if(r1!=null){
								int x=TrimRead.trimFast(r1, trimLeft, trimRight, trimq, 1);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2!=null){
								int x=TrimRead.trimFast(r2, trimLeft, trimRight, trimq, 1);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
						}
					}

					if(minReadLength>0 && !remove){
						int rlen=(r1==null || r1.bases==null ? 0 : r1.bases.length);
						int rlen2=(r2==null || r2.bases==null ? 0 : r2.bases.length);
						if(rlen<minReadLength && rlen2<minReadLength){
							basesTrimmedT+=(rlen+rlen2);
							remove=true;
						}
					}

					if(!remove){
						if(r2!=null){r2.reverseComplement();}
						readsProcessed++;

						final int[] rvector=new int[5];
						int trueSize=r1.insertSizeMapped(ignoreMappingStrand);
						if(verbose){System.err.println("True Insert: "+trueSize);}

						int bestInsert;
						int bestScore=-1;
						int bestGood=-1;
						int bestBad=999999;
						boolean ambig, tooShort=false;

						//					assert(false) : r+"\n"+(USE_MAPPING)+", "+(r.chrom==r.mate.chrom)+", "+()+", "+()+", "+()+", "+()+", ";

						if(r2==null){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=r1.bases.length;
							assert(r1.bases.length==r1.insert()) : r1.bases.length+" != "+r1.insert()+"; actual = "+trueSize;
							//						if(bestInsert!=trueSize){
							//							System.err.println("Bad insert size for pre-joined read "+r.numericID+": len="+r.bases.length+", insert="+r.insert()+", actual="+trueSize);
							//						}
							ambig=false;
						}else if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && ((r1.mapped() || r1.synthetic()) && (r2.mapped() || r2.synthetic()))){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=trueSize;
							ambig=false;
						}else if(SKIP_MATED_READS && r1.insertvalid() && r1.insert()>0){
							bestScore=100;
							bestGood=30;
							bestBad=0;
							bestInsert=r1.insert();
							ambig=false;
						}else{
							if(MATE_BY_OVERLAP){
								bestInsert=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0, MIN_OVERLAPPING_BASES);
								bestScore=rvector[0];
								bestGood=rvector[1];
								bestBad=rvector[2];
								ambig=(rvector[4]==1);
								final int len1=r1.bases.length, len2=r2.bases.length;
								for(int trims=0, q=trimq; trims<TRIM_ON_OVERLAP_FAILURE && !qtrim && bestInsert<0 /*&& !ambig*/; trims++, q+=8){
//									System.err.println(trims+", "+q);
									Object old1=r1.obj;
									Object old2=r2.obj;
									tr1=TrimRead.trim(r1, false, true, q, 1+len1*4/10); //r1.bases.length);
									tr2=TrimRead.trim(r2, true, false, q, 1+len2*4/10); //r2.bases.length);
									r1.obj=old1;
									r2.obj=old2;
									if(tr1!=null || tr2!=null){
//										System.err.println(r1.bases.length+", "+r2.bases.length);
										int x=BBMergeOverlapper.mateByOverlap(r1, r2, rvector, MIN_OVERLAPPING_BASES_0-1, MIN_OVERLAPPING_BASES);
										if(x>-1){
//											System.err.println(trims);
											bestInsert=x;
											bestScore=rvector[0];
											bestGood=rvector[1];
											bestBad=rvector[2];
											ambig=(rvector[4]==1);
											trims=TRIM_ON_OVERLAP_FAILURE;
										}else{
											if(tr1!=null){tr1.untrim();}
											if(tr2!=null){tr2.untrim();}
										}
									}
								}
							}else{
								ambig=false;
								bestInsert=-1;
							}
						}

						tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
						
						if(joinperfectonly && bestBad>0){ambig=true;}

						if(ambig){ambiguousCount++;}
						else if(tooShort){tooShortCount++;}
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
							
							if((ambig || bestInsert<0 || tooShort) && (rosbad!=null || !MIX_BAD_AND_GOOD)){
								if(rosbad!=null){
									if(listb!=null){listb.add(r1);}
								}
							}else{
								if(listg!=null){
									Read x=r1;
									if(joinReads && r2!=null){x=r1.joinRead(bestInsert);}
									listg.add(x);
								}
							}

							if(rosinsert!=null){
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
									sb.append("\t"+bestGood+"\t"+bestBad+"\t"+bestScore);
								}
								r1.obj=sb;
								listi.add(r1);
							}
						}
						if(r2!=null){r2.reverseComplement();}
					}
				}

				if(rosgood!=null){rosgood.add(listg, ln.id);}
				if(rosbad!=null){rosbad.add(listb, ln.id);}
				if(rosinsert!=null){rosinsert.add(listi, ln.id);}

				//			System.err.println("returning list");
				cris.returnList(ln, ln.list.isEmpty());
				//			System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			System.err.println("reads: "+(reads==null ? "null" : reads.size()));
			}
			cris.returnList(ln, ln.list.isEmpty());
		}

		int[] hist=new int[1000];

		long readsProcessed=0;
		long matedCount=0;
		long correctCount=0;
		long ambiguousCount=0;
		long tooShortCount=0;
		long incorrectCount=0;
		long noSolutionCount=0;
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		int insertMax=0;
		int insertMin=999999999;
		
		long basesTrimmedT=0;
		long readsTrimmedT=0;
		
		private final ConcurrentReadStreamInterface cris;
		private final RTextOutputStream3 rosgood;
		private final RTextOutputStream3 rosbad;
		private final RTextOutputStream3 rosinsert;
		
		private final boolean joinReads;
		private final boolean joinperfectonly;
	}
	
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
	
	private long maxReads=-1;
	private boolean join=true;
	
	private byte qin=-1;
	private byte qout=-1;
	
	static boolean errorState=false;
	
	static boolean trimRight=false;
	static boolean trimLeft=false;
	static boolean untrim=false;
	static byte trimq=6;
	static int minReadLength=0;
	static int minInsert=0;
	static boolean qtrim=false;
	static int TRIM_ON_OVERLAP_FAILURE=1;
	
	static int[] histTotal=new int[1000];
	static int bin=1;

	static long readsProcessedTotal=0;
	static long matedCountTotal=0;
	static long correctCountTotal=0;
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long incorrectCountTotal=0;
	static long noSolutionCountTotal=0;
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	static long basesTrimmedTotal=0;
	static long readsTrimmedTotal=0;
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;
	
	private static int MIN_OVERLAPPING_BASES=12;
	private static int MIN_OVERLAPPING_BASES_0=8;
	
	private static boolean MATE_BY_OVERLAP=true;
	private static boolean SKIP_MATED_READS=false;
	private static boolean OUTPUT_FAILED=true;
	private static boolean MIX_BAD_AND_GOOD=false;
	
	private static boolean overwrite=true;
	private static boolean append=false;
	private static boolean verbose=false;
	
	private static int THREADS=-1;
	private static float version=3.0f;
	
}
