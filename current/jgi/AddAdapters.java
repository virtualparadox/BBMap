package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;
import align2.ListNum;
import align2.QualityTools;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;

/**
 * @author Brian Bushnell
 * @date Mar 16, 2014
 *
 */
public class AddAdapters {
	


	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		AddAdapters rr=new AddAdapters(args);
		if(rr.writeMode){
			rr.write(t);
		}else{
			rr.read(t);
		}
	}
	
	private void printOptions(){
		System.err.println("Please consult the shellscript for usage information.");
	}
	
	public AddAdapters(String[] args){
		
		if(args==null || args.length==0){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");

		boolean setxs=false, setintron=false;
		boolean setInterleaved=false; //Whether it was explicitly set.

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.READ_BUFFER_NUM_BUFFERS=Tools.min(8, Shared.READ_BUFFER_NUM_BUFFERS);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("null") || a.equals(in2)){
				// do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.THREADS=Tools.max(Integer.parseInt(b), 1);
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("adapter") || a.equals("adapters") || a.equals("ref")){
				adapterFile=b;
			}else if(a.equals("literal") || a.equals("literals")){
				literals=(b==null ? null : b.split(","));
			}else if(a.equals("minlength") || a.equals("minlen") || a.equals("ml")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("3'") || a.equalsIgnoreCase("3prime") || a.equalsIgnoreCase("3-prime") || a.equalsIgnoreCase("right") || a.equalsIgnoreCase("r")){
				right=Tools.parseBoolean(b);
			}else if(a.equals("5'") || a.equalsIgnoreCase("5prime") || a.equalsIgnoreCase("5-prime") || a.equalsIgnoreCase("left") || a.equalsIgnoreCase("l")){
				right=!Tools.parseBoolean(b);
			}else if(a.equals("end")){
				if(b.equals("3'") || b.equalsIgnoreCase("3prime") || b.equalsIgnoreCase("3-prime") || b.equalsIgnoreCase("right") || a.equalsIgnoreCase("r")){
					right=true;
				}if(b.equals("5'") || b.equalsIgnoreCase("5prime") || b.equalsIgnoreCase("5-prime") || b.equalsIgnoreCase("left") || a.equalsIgnoreCase("l")){
					right=true;
				}
			}else if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription")){
				Shared.TRIM_READ_COMMENTS=Tools.parseBoolean(b);
			}else if(a.equals("addslash")){
				addslash=Tools.parseBoolean(b);
			}else if(a.equals("adderrors")){
				adderrors=Tools.parseBoolean(b);
			}else if(a.equals("addreversecomplement") || a.equals("arc")){
				addRC=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("write")){
				writeMode=Tools.parseBoolean(b);
			}else if(a.equals("grade")){
				writeMode=!Tools.parseBoolean(b);
			}else if(a.equals("mode")){
				if("grade".equalsIgnoreCase(b) || "read".equalsIgnoreCase(b)){
					writeMode=false;
				}else if("generate".equalsIgnoreCase(b) || "write".equalsIgnoreCase(b) || "add".equalsIgnoreCase(b)){
					writeMode=true;
				}else{
					throw new RuntimeException("Unknown mode "+b);
				}
			}else if(a.equals("fastareadlen") || a.equals("fastareadlength")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.equals("fastaminread") || a.equals("fastaminlen") || a.equals("fastaminlength")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
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
			}else if(a.equals("qauto")){
				FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
			}else if(a.equals("tuc") || a.equals("touppercase")){
				Read.TO_UPPER_CASE=Tools.parseBoolean(b);
			}else if(a.equals("lctn") || a.equals("lowercaseton")){
				Read.LOWER_CASE_TO_N=Tools.parseBoolean(b);
			}else if(a.equals("tossbrokenreads") || a.equals("tbr")){
				boolean x=Tools.parseBoolean(b);
				Read.NULLIFY_BROKEN_QUALITY=x;
				ConcurrentGenericReadInputStream.REMOVE_DISCARDED_READS=x;
			}else if(a.equals("testinterleaved")){
				FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				outstream.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
				setInterleaved=true;
			}else if(a.equals("forceinterleaved")){
				FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
				outstream.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				setInterleaved=true;
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
					setInterleaved=true;
				}
			}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(out1==null && i==1 && !arg.contains("=")){
				out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
//		if(maxReads!=-1){ReadWrite.USE_GUNZIP=ReadWrite.USE_UNPIGZ=false;}
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.THREADS>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(writeMode && out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
			out1="stdout.fq";
		}
		
		if(!setInterleaved){
			assert(in1!=null && (!writeMode || out1!=null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else if(writeMode){ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, false, out1, out2)){
			throw new RuntimeException("\n\nOVERWRITE="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		ffa=FileFormat.testInput(adapterFile, FileFormat.FASTA, null, true, true);
		
		adapters=makeAdapterList();
		
		if(writeMode){
			if(adapters==null || adapters.isEmpty()){
				throw new RuntimeException("\n\nPlease specify adapters with 'adapters=file.fa' or 'literal=AGCTACGT'\n");
			}
			randy=new Random();
		}
	}
	
	private final ArrayList<byte[]> makeAdapterList(){
		if(ffa==null && literals==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		if(ffa!=null){
			FastaReadInputStream fris=new FastaReadInputStream(ffa, false, false, -1);
			for(Read r=fris.next(); r!=null; r=fris.next()){
				if(r.bases!=null){
					list.add(r.bases);
				}
			}
			fris.close();
		}
		if(literals!=null){
			for(String s : literals){
				if(s!=null && !"null".equalsIgnoreCase(s)){
					list.add(s.getBytes());
				}
			}
		}
		
		if(addRC){
			int x=list.size();
			for(int i=0; i<x; i++){
				list.add(AminoAcid.reverseComplementBases(list.get(i)));
			}
		}
		
		return list.size()>0 ? list : null;
	}
	
	void write(Timer t){
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, false, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		RTextOutputStream3 ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=new RTextOutputStream3(ffout1, ffout2, null, null, buff, null, false);
			ros.start();
		}

		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=(r1.bases==null ? 0 : r1.bases.length);
					final int initialLength2=(r2==null ? 0 : r2.bases==null ? 0 : r2.bases.length);
					
					{
						addAdapter(r1);
					}
					if(r2!=null){
						addAdapter(r2);
					}
					
					if(r2==null){
						r1.id=r1.numericID+"_"+r1.id;
					}else{
						String base=r1.numericID+"_"+r1.id+"_"+r2.id;
						if(addslash){
							r1.id=base+" /1";
							r2.id=base+" /2";
						}else{
							r1.id=base;
							r2.id=base;
						}
					}
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
//		System.err.println(cris.errorState()+", "+(ros==null ? "null" : (ros.errorState()+", "+ros.finishedSuccessfully())));
//		if(ros!=null){
//			ReadStreamWriter rs1=ros.getRS1();
//			ReadStreamWriter rs2=ros.getRS2();
//			System.err.println(rs1==null ? "null" : rs1.finishedSuccessfully());
//			System.err.println(rs2==null ? "null" : rs2.finishedSuccessfully());
//		}
//		assert(false);
		
		t.stop();

		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		outstream.println("Adapters Added:         \t"+adaptersAdded+" reads ("+String.format("%.2f",adaptersAdded*100.0/readsProcessed)+"%) \t"+
				adapterBasesAdded+" bases ("+String.format("%.2f",adapterBasesAdded*100.0/basesProcessed)+"%)");

		outstream.println("Valid Output:           \t"+validReads+" reads ("+String.format("%.2f",validReads*100.0/readsProcessed)+"%) \t"+
				validBases+" bases ("+String.format("%.2f",validBases*100.0/basesProcessed)+"%)");

		
		outstream.println("\nTime:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void addAdapter(Read r){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final int remaining, initial=(bases==null ? 0 : bases.length);
		final byte[] adapter;
		final int loc;
		int ab=0, rb=0;
		
		readsProcessed++;
		basesProcessed+=initial;
		if(initial>0 && randy.nextFloat()<adapterProb){
			adapter=adapters.get(randy.nextInt(adapters.size()));
			loc=randy.nextInt(initial);
			adaptersAdded++;

			if(right){
				final int lim=Tools.min(initial, adapter.length+loc);
				for(int i=loc, j=0; i<lim; i++, j++){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=adapter[j];
						if(adderrors){
							byte q=(quals==null ? 30 : quals[i]);
							if(randy.nextFloat()<QualityTools.PROB_ERROR[q]){
								int old=AminoAcid.baseToNumber[bases[i]];
								bases[i]=AminoAcid.numberToBase[(old+randy.nextInt(3))%4];
							}
						}
					}
					ab++;
				}
				for(int i=lim; i<initial; i++){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=AminoAcid.numberToBase[randy.nextInt(4)];
					}
					rb++;
				}
				remaining=loc;
			}else{
				final int lim=Tools.max(-1, loc-adapter.length);
				for(int i=loc, j=adapter.length-1; i>lim; i--, j--){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=adapter[j];
						if(adderrors){
							byte q=(quals==null ? 30 : quals[i]);
							if(randy.nextFloat()<QualityTools.PROB_ERROR[q]){
								int old=AminoAcid.baseToNumber[bases[i]];
								bases[i]=AminoAcid.numberToBase[(old+randy.nextInt(3))%4];
							}
						}
					}
					ab++;
				}
				for(int i=lim; i>-1; i--){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=AminoAcid.numberToBase[randy.nextInt(4)];
					}
					rb++;
				}
				remaining=initial-loc-1;
			}
			assert(remaining<initial) : "\nremaining="+remaining+", initial="+initial+", rb="+rb+", ab="+ab+
				", loc="+loc+", adapter.length="+(adapter==null ? 0 : adapter.length)+"\n";
		}else{
			adapter=null;
			loc=-1;
			remaining=initial;
		}
		
		assert(remaining==initial-(rb+ab));
		assert(remaining>=0);

		adapterBasesAdded+=ab;
		randomBasesAdded+=rb;
		r.id=initial+"_"+remaining;
		if(remaining>=minlen){
			validReads++;
			validBases+=remaining;
		}
	}
	
	/*--------------------------------------------------------------*/
	
	void read(Timer t){

		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, false, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					grade(r1, r2);
				}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStream(cris);
		
		t.stop();
		
		long validBasesRemoved=validBasesExpected-validBasesCounted;
		long incorrect=readsProcessed-correct;
		long incorrectBases=basesProcessed-correctBases;

		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Total output:                        \t"+readsProcessed+" reads                  \t"+basesProcessed+" bases          ");
		outstream.println("Perfectly Correct (% of output):     \t"+correct+" reads ("+String.format("%.3f",correct*100.0/readsProcessed)+
				"%)        \t"+correctBases+" bases ("+String.format("%.3f",correctBases*100.0/basesProcessed)+"%)");
		outstream.println("Incorrect (% of output):             \t"+incorrect+" reads ("+String.format("%.3f",incorrect*100.0/readsProcessed)+
				"%)        \t"+incorrectBases+" bases ("+String.format("%.3f",incorrectBases*100.0/basesProcessed)+"%)");
		outstream.println();
//		outstream.println("Too Short:              \t"+tooShort+" reads ("+String.format("%.3f",tooShort*100.0/readsProcessed)+"%) \t"+
//				tooShortBases+" bases ("+String.format("%.3f",tooShortBases*100.0/basesProcessed)+"%)");
//		outstream.println("Too Long:               \t"+tooLong+" reads ("+String.format("%.3f",tooLong*100.0/readsProcessed)+"%) \t"+
//				tooLongBases+" bases ("+String.format("%.3f",tooLongBases*100.0/basesProcessed)+"%)");
		
		outstream.println("Adapters Remaining (% of adapters):  \t"+(adapterReadsRemaining)+" reads ("+String.format("%.3f",adapterReadsRemaining*100.0/adapterReadsTotal)+
				"%)        \t"+adapterBasesRemaining+" bases ("+String.format("%.3f",adapterBasesRemaining*100.0/basesProcessed)+"%)");
		outstream.println("Non-Adapter Removed (% of valid):    \t"+tooShort+" reads ("+String.format("%.3f",tooShort*100.0/readsProcessed)+
				"%)        \t"+validBasesRemoved+" bases ("+String.format("%.3f",validBasesRemoved*100.0/validBasesExpected)+"%)");
		
		if(broken>0 || mispaired>0){
			outstream.println("Broken:                 \t"+broken+" reads ("+String.format("%.2f",broken*100.0/readsProcessed)+"%");
			outstream.println("Mispaired:              \t"+mispaired+" reads ("+String.format("%.2f",mispaired*100.0/readsProcessed)+"%");
		}

		
//		outstream.println("\nTime:                         \t"+t);
//		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
//		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void grade(Read r1, Read r2){
		final String a=r1.id.split(" ")[0];
		final String b=(r2==null ? a : r2.id.split(" ")[0]);
		final int len=a.split("_").length;
		
		if(r2!=null){
			if(len!=5){
				throw new RuntimeException("Headers are corrupt, or unpaired reads are being processed as paired.  Try running with 'int=f'");
			}
			if(r1.id.endsWith(" /2") || r2.id.endsWith(" /1") || !a.equals(b)){
				mispaired+=2;
			}
			if(r1.id.endsWith(" /2")){r1.setPairnum(1);}
			if(r2.id.endsWith(" /1")){r2.setPairnum(0);}
		}else{
			if(len!=3){
				throw new RuntimeException("Headers are corrupt, or paired reads are being processed as unpaired.  Try running with 'int=t' or with 'in1=' and 'in2='.");
			}
		}
		grade(r1);
		grade(r2);
	}
	
	private void grade(Read r){
		if(r==null){return;}
		final int offset=(2*r.pairnum());
		
		String[] sa=r.id.split(" ")[0].split("_");
		final long id=Long.parseLong(sa[0]);
		final int initial=Integer.parseInt(sa[1+offset]);
		final int remaining=Integer.parseInt(sa[2+offset]);
		final int actual=r.bases==null ? 0 : r.bases.length;
		
		readsProcessed++;
		basesProcessed+=actual;
		
		assert(initial>=remaining);
		
		if(actual>initial){broken++;}
		
		validBasesExpected+=remaining;
		
		if(initial==remaining){//Should not have trimmed
			if(actual==remaining || (actual<2 && (remaining<1 || remaining<minlen))){
				correct++;
				correctBases+=remaining;
				validBasesCounted+=remaining;
				trueNeg++;
			}else if(actual<remaining){
				tooShort++;
				tooShortReadBases+=actual;
				tooShortBases+=(remaining-actual);
				validBasesCounted+=actual;
				falsePos++;
			}else if(actual>remaining){
				tooLong++;
				tooLongReadBases+=remaining;
				tooLongBases+=(actual-remaining);
				validBasesCounted+=remaining;
				falseNeg++;
			}
		}else{//Should have trimmed
			
			adapterBasesTotal+=(initial-remaining);
			adapterReadsTotal++;
			
			if(actual==remaining || (actual<2 && (remaining<1 || remaining<minlen))){
				correct++;
				correctBases+=remaining;
				validBasesCounted+=remaining;
				truePos++;
			}else if(actual<remaining){
				tooShort++;
				tooShortReadBases+=actual;
				tooShortBases+=(remaining-actual);
				validBasesCounted+=actual;
				truePos++;
			}else if(actual>remaining){
				tooLong++;
				tooLongReadBases+=actual;
				tooLongBases+=(actual-remaining);
				adapterBasesRemaining+=(actual-remaining);
				validBasesCounted+=remaining;
				falseNeg++;
				adapterReadsRemaining++;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	public boolean errorState=false;
	
	private String in1=null;
	private String in2=null;

	private String out1=null;
	private String out2=null;
	
	private String extin=null;
	private String extout=null;

	private String adapterFile=null;
	private String[] literals=null;
	
	private boolean overwrite=false;

	/** Add /1 and /2 to paired reads */
	private boolean addslash=true;
	/** Encode correct answer in read ID field */
	private boolean changename=true;
	/** Add errors from quality value */
	private boolean adderrors=true;

	private boolean addRC=true;
	/** aka 3' */
	private boolean right=true;
	
	private long maxReads=-1;
	private int minlen=1;
	
	private byte qin=-1;
	private byte qout=-1;
	
	private boolean writeMode=true;
	private float adapterProb=0.5f;
	
	private long readsProcessed=0;
	private long basesProcessed=0;
	private long adaptersAdded=0;
	private long adapterBasesAdded=0;
	private long randomBasesAdded=0;
	private long validReads=0;
	private long validBases=0;

	private long truePos=0;
	private long trueNeg=0;
	private long falsePos=0;
	private long falseNeg=0;
	private long broken=0;
	private long mispaired=0;
	
	private long tooShort=0;
	private long tooLong=0;
	private long correct=0;
	private long fullyRemoved=0;

	private long tooShortBases=0;
	private long tooLongBases=0;
	private long tooShortReadBases=0;
	private long tooLongReadBases=0;
	private long correctBases=0;

	private long validBasesCounted=0;
	private long validBasesExpected=0;
	
//	private long invalidBasesCounted=0;
	private long adapterBasesTotal=0;
	private long adapterReadsTotal=0;
	private long adapterReadsRemaining=0;
	private long adapterBasesRemaining=0;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	private final FileFormat ffa;
	
	private final ArrayList<byte[]> adapters;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	
	private java.util.Random randy;
	
}
