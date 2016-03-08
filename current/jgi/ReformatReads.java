package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;
import stream.SamLine;

import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;

import align2.ListNum;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;

/**
 * @author Brian Bushnell
 * @date Sep 11, 2012
 *
 */
public class ReformatReads {

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		ReformatReads rr=new ReformatReads(args);
		rr.process(t);
	}
	
	public ReformatReads(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.READ_BUFFER_NUM_BUFFERS=Tools.min(8, Shared.READ_BUFFER_NUM_BUFFERS);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("null") || a.equals(parser.in2)){
				// do nothing
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("sample") || a.equals("samplereads") || a.equals("samplereadstarget") || a.equals("srt")){
				sampleReadsTarget=Long.parseLong(b);
				sampleReadsExact=(sampleReadsTarget>0);
			}else if(a.equals("samplebases") || a.equals("samplebasestarget") || a.equals("sbt")){
				sampleBasesTarget=Long.parseLong(b);
				sampleBasesExact=(sampleBasesTarget>0);
			}else if(a.equals("parsecustom")){
				parsecustom=Tools.parseBoolean(b);
			}else if(a.equals("addslash")){
				addslash=Tools.parseBoolean(b);
			}else if(a.equals("uniquenames")){
				uniqueNames=Tools.parseBoolean(b);
			}else if(a.equals("verifyinterleaved") || a.equals("verifyinterleaving") || a.equals("vint")){
				verifyinterleaving=Tools.parseBoolean(b);
			}else if(a.equals("rcompmate") || a.equals("rcm")){
				reverseComplimentMate=Tools.parseBoolean(b);
				outstream.println("Set RCOMPMATE to "+reverseComplimentMate);
			}else if(a.equals("rcomp") || a.equals("rc")){
				reverseCompliment=Tools.parseBoolean(b);
				outstream.println("Set RCOMP to "+reverseCompliment);
			}else if(a.equals("deleteempty") || a.equals("deletempty") || a.equals("delempty") || a.equals("def")){
				deleteEmptyFiles=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					parser.in1=b.replace("#", "1");
					parser.in2=b.replace("#", "2");
				}
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Download parser fields
			
			maxReads=parser.maxReads;	
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;

			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;

			trimq=parser.trimq;
			minAvgQuality=parser.minAvgQuality;
			averageQualityByProbability=parser.averageQualityByProbability;
			minReadLength=parser.minReadLength;
			maxReadLength=parser.maxReadLength;
			minLenFraction=parser.minLenFraction;
			breakLength=parser.breakLength;
			requireBothBad=parser.requireBothBad;
			trimBadSequence=parser.trimBadSequence;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			testsize=parser.testsize;
			

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			outsingle=parser.outsingle;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		if(SamLine.setxs && !SamLine.setintron){SamLine.INTRON_LIMIT=10;}
		qtrim=qtrimLeft||qtrimRight;

		if(verifyinterleaving){
			setInterleaved=true;
//			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
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
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.THREADS>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
			if(!parser.setOut){
				out1="stdout";
			}
		}
		
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		if(outsingle!=null && outsingle.equalsIgnoreCase("null")){outsingle=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outsingle)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		FASTQ.PARSE_CUSTOM=parsecustom;
		
		{
			byte qin=Parser.qin, qout=Parser.qout;
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
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutsingle=FileFormat.testOutput(outsingle, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

//		System.err.println("\n"+ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1+"\n");
//		assert(false) : ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1;

		if(ffin1!=null && ffout1!=null && ffin1.samOrBam()){
			if(ffout1.samOrBam()){
				useSharedHeader=true;
				SamLine.CONVERT_CIGAR_TO_MATCH=true;
			}else if(ffout1.bread()){
				SamLine.CONVERT_CIGAR_TO_MATCH=true;
			}
		}

		nameMap1=(uniqueNames ? new HashMap<String, Integer>() : null);
		nameMap2=(uniqueNames ? new HashMap<String, Integer>() : null);
	}
	
	void process(Timer t){
		
		long readsRemaining=0;
		long basesRemaining=0;
		
		if(sampleReadsExact || sampleBasesExact){
			long[] counts=countReads(in1, in2, maxReads);
			readsRemaining=counts[0];
			basesRemaining=counts[2];
			setSampleSeed(sampleseed);
		}
		
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
//		if(verbose){
			System.err.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
//		}
		
		assert(!paired || breakLength<1) : "Paired input cannot be broken with 'breaklength'";

		final RTextOutputStream3 ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=new RTextOutputStream3(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader);
			ros.start();
		}else{ros=null;}
		
		final RTextOutputStream3 rosb;
		if(outsingle!=null){
			final int buff=4;
			
			rosb=new RTextOutputStream3(ffoutsingle, null, buff, null, useSharedHeader);
			rosb.start();
		}else{rosb=null;}
		final boolean discardTogether=(!paired || (outsingle==null && !requireBothBad));
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		//Only used with deleteEmptyFiles flag
		long readsOut1=0;
		long readsOut2=0;
		long readsOutSingle=0;
		
		long basesTrimmedT=0;
		long readsTrimmedT=0;
		
		long lowqBasesT=0;
		long lowqReadsT=0;
		
		long readShortDiscardsT=0;
		long baseShortDiscardsT=0;

		final boolean MAKE_QHIST=ReadStats.COLLECT_QUALITY_STATS;
		final boolean MAKE_QAHIST=ReadStats.COLLECT_QUALITY_ACCURACY;
		final boolean MAKE_MHIST=ReadStats.COLLECT_MATCH_STATS;
		final boolean MAKE_BHIST=ReadStats.COLLECT_BASE_STATS;
		final ReadStats readstats=(MAKE_QHIST || MAKE_MHIST || MAKE_BHIST) ? new ReadStats() : null;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){
				ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(32));
				
				if(breakLength>0){
					breakReads(reads, breakLength, minReadLength);
				}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=(r1.bases==null ? 0 : r1.bases.length);
					final int initialLength2=(r2==null ? 0 : r2.bases==null ? 0 : r2.bases.length);

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					if(readstats!=null){
						if(MAKE_QHIST){readstats.addToQualityHistogram(r1);}
						if(MAKE_BHIST){readstats.addToBaseHistogram(r1);}
						if(MAKE_MHIST){readstats.addToMatchHistogram(r1);}
						if(MAKE_QAHIST){readstats.addToQualityAccuracy(r1);}
					}
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
						if(reverseCompliment){r1.reverseComplement();}
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
						if(reverseCompliment || reverseComplimentMate){r2.reverseComplement();}
					}
					
					if(verifyinterleaving){
						String s1=r1==null ? null : r1.id;
						String s2=r2==null ? null : r2.id;
						boolean b=FASTQ.testPairNames(s1, s2);
						if(!b){
							outstream.println("Names do not appear to be correctly paired.\n"+s1+"\n"+s2+"\n");
							ReadWrite.closeStreams(cris, ros);
							System.exit(1);
						}
					}
					
					if(trimBadSequence){//Experimental
						if(r1!=null){
							int x=TrimRead.trimBadSequence(r1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null){
							int x=TrimRead.trimBadSequence(r2);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
					}

					if(minAvgQuality>0){
						final int a=(r1==null ? 99 : averageQualityByProbability ? r1.avgQualityByProbability() : r1.avgQuality());
						final int b=(r2==null ? 99 : averageQualityByProbability ? r2.avgQualityByProbability() : r2.avgQuality());
						if((a<minAvgQuality || b<minAvgQuality)){
							if(r1!=null && a<minAvgQuality && !r1.discarded()){
								lowqBasesT+=r1.bases.length;
								lowqReadsT++;
								r1.setDiscarded(true);
							}
							if(r2!=null && b<minAvgQuality && !r2.discarded()){
								lowqBasesT+=r2.bases.length;
								lowqReadsT++;
							}
						}
					}
					
					if(forceTrimLeft>0 || forceTrimRight>0){
						if(r1!=null && !r1.discarded()){
							int x=TrimRead.trimToPosition(r1, forceTrimLeft>0 ? forceTrimLeft : 0, forceTrimRight>0 ? forceTrimRight : r1.bases.length, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null && !r2.discarded()){
							int x=TrimRead.trimToPosition(r2, forceTrimLeft>0 ? forceTrimLeft : 0, forceTrimRight>0 ? forceTrimRight : r2.bases.length, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(qtrim){
						if(r1!=null && !r1.discarded()){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null && !r2.discarded()){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(minlen1>0 || minlen2>0 || maxReadLength>0){
//						assert(false) : minlen1+", "+minlen2+", "+maxReadLength+", "+r1.bases.length;
						if(r1!=null && !r1.discarded()){
							int rlen=(r1.bases==null ? 0 : r1.bases.length);
							if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
								r1.setDiscarded(true);
								readShortDiscardsT++;
								baseShortDiscardsT+=rlen;
							}
						}
						if(r2!=null && !r2.discarded()){
							int rlen=(r2.bases==null ? 0 : r2.bases.length);
							if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
								r2.setDiscarded(true);
								readShortDiscardsT++;
								baseShortDiscardsT+=rlen;
							}
						}
					}
					
					boolean remove=false;
					if(r2==null){
						remove=r1.discarded();
					}else{
						remove=requireBothBad ? (r1.discarded() && r2.discarded()) : (r1.discarded() || r2.discarded());
					}
					
					if(remove){reads.set(idx, null);}
					else{
						if(uniqueNames){
							{
								Integer v=nameMap1.get(r1.id);
								if(v==null){
									nameMap1.put(r1.id, 1);
								}else{
									v++;
									nameMap1.put(r1.id, v);
									r1.id=r1.id+"_"+v;
								}
							}
							if(r2!=null){
								Integer v=nameMap2.get(r2.id);
								if(v==null){
									nameMap2.put(r2.id, 1);
								}else{
									v++;
									nameMap2.put(r2.id, v);
									r2.id=r2.id+"_"+v;
								}
							}
						}
						if(addslash){
							if(r1.id==null){r1.id=" "+r1.numericID;}
							if(!r1.id.contains(" /1")){r1.id+=" /1";}
							if(r2!=null){
								if(r2.id==null){r2.id=" "+r2.numericID;}
								if(!r2.id.contains(" /2")){r2.id+=" /2";}
							}
						}
					}
					
					if(singles!=null){
						if(r1.discarded() || (r2!=null && r2.discarded())){
							if(!r1.discarded()){
								Read r=r1.clone();
								r.mate=null;
								r.setPairnum(0);
								singles.add(r);
							}else if(r2!=null && !r2.discarded()){
								Read r=r2.clone();
								r.mate=null;
								r.setPairnum(0);
								singles.add(r);
							}
						}
					}
				}
				
				final ArrayList<Read> listOut;
				
//				assert(false) : sampleReadsExact+", "+sampleBasesExact;
				if(sampleReadsExact || sampleBasesExact){
					listOut=new ArrayList<Read>();
					if(sampleReadsExact){
						for(Read r : reads){
							if(r!=null){
								assert(readsRemaining>0) : readsRemaining;
								double prob=sampleReadsTarget/(double)(readsRemaining);
//								System.err.println("sampleReadsTarget="+sampleReadsTarget+", readsRemaining="+readsRemaining+", prob="+prob);
								if(randy.nextDouble()<prob){
									listOut.add(r);
									sampleReadsTarget--;
								}
							}
							readsRemaining--;
						}
					}else if(sampleBasesExact){
						for(Read r : reads){
							if(r!=null){
								assert(basesRemaining>0) : basesRemaining;
								int bases=r.bases.length+(r.mate==null ? 0 : r.mate.bases.length);
								double prob=sampleBasesTarget/(double)(basesRemaining);
								if(randy.nextDouble()<prob){
									listOut.add(r);
									sampleBasesTarget-=bases;
								}
								basesRemaining-=bases;
							}
						}
					}
				}else{
					listOut=reads;
				}
				if(deleteEmptyFiles){
					for(Read r : listOut){
						if(r!=null){
							readsOut1++;
							if(r.mate!=null){
								readsOut2++;
							}
						}
					}
					if(singles!=null){readsOutSingle+=singles.size();}
				}
				if(ros!=null){ros.add(listOut, ln.id);}
				if(rosb!=null){rosb.add(singles, ln.id);}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadStats.writeAll(paired);
		
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
		
		if(deleteEmptyFiles){
			deleteEmpty(readsOut1, readsOut2, readsOutSingle);
		}
		
//		System.err.println(cris.errorState()+", "+(ros==null ? "null" : (ros.errorState()+", "+ros.finishedSuccessfully())));
//		if(ros!=null){
//			ReadStreamWriter rs1=ros.getRS1();
//			ReadStreamWriter rs2=ros.getRS2();
//			System.err.println(rs1==null ? "null" : rs1.finishedSuccessfully());
//			System.err.println(rs2==null ? "null" : rs2.finishedSuccessfully());
//		}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		if(qtrim || forceTrimLeft>0 || forceTrimRight>0 || trimBadSequence){
			outstream.println("QTrimmed:               \t"+readsTrimmedT+" reads ("+String.format("%.2f",readsTrimmedT*100.0/readsProcessed)+"%) \t"+
					basesTrimmedT+" bases ("+String.format("%.2f",basesTrimmedT*100.0/basesProcessed)+"%)");
		}else if(minReadLength>0 || maxReadLength>0){
			outstream.println("Short Read Discards:    \t"+readShortDiscardsT+" reads ("+String.format("%.2f",readShortDiscardsT*100.0/readsProcessed)+"%) \t"+
					baseShortDiscardsT+" bases ("+String.format("%.2f",baseShortDiscardsT*100.0/basesProcessed)+"%)");
		}
		if(minAvgQuality>0){
			outstream.println("Low quality discards:   \t"+lowqReadsT+" reads ("+String.format("%.2f",lowqReadsT*100.0/readsProcessed)+"%) \t"+
					lowqBasesT+" bases ("+String.format("%.2f",lowqBasesT*100.0/basesProcessed)+"%)");
		}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length())+
					(qfin1==null ? 0 : new File(qfin1).length())+(qfin2==null ? 0 : new File(qfin2).length()));//*passes
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+String.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(verifyinterleaving){
			outstream.println("Names appear to be correctly paired.");
		}
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void deleteEmpty(long readsOut1, long readsOut2, long readsOutSingle){
		deleteEmpty(readsOut1, ffout1, qfout1);
		deleteEmpty(readsOut2, ffout2, qfout2);
		deleteEmpty(readsOutSingle, ffoutsingle, null);
	}
	
	private void deleteEmpty(long count, FileFormat ff, String qf){
		try {
			if(ff!=null && count<1){
				String s=ff.name();
				if(s!=null && !ff.stdio() && !ff.devnull()){
					File f=new File(ff.name());
					if(f.exists()){
						f.delete();
					}
				}
				if(qf!=null){
					File f=new File(qf);
					if(f.exists()){
						f.delete();
					}
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

	public static void breakReads(ArrayList<Read> list, final int max, int min){
		if(!containsReadsAboveSize(list, max)){return;}
		assert(max>0) : "max read length must be positive.";
		assert(max>=min) : "max read length must be at least min read length: "+max+"<"+min;
		min=Tools.max(0, min);
		
		ArrayList<Read> temp=new ArrayList<Read>(list.size()*2);
		for(Read r : list){
			if(r==null || r.bases==null || r.bases.length<=max){
				temp.add(r);
			}else{
				final byte[] bases=r.bases;
				final byte[] quals=r.quality;
				final String name=r.id;
				final int limit=bases.length-min;
				for(int num=1, start=0, stop=max; start<limit; num++, start+=max, stop+=max){
					if(verbose){
						System.err.println(bases.length+", "+start+", "+stop);
						if(quals!=null){System.err.println(quals.length+", "+start+", "+stop);}
					}
					stop=Tools.min(stop, bases.length);
					byte[] b2=Arrays.copyOfRange(bases, start, stop);
					byte[] q2=(quals==null ? null : Arrays.copyOfRange(quals, start, stop));
					String n2=name+"_"+num;
					Read r2=new Read(b2, -1, -1, -1, n2, q2, r.numericID, r.flags);
					r2.setMapped(false);
					temp.add(r2);
				}
			}
		}
		list.clear();
		list.ensureCapacity(temp.size());
		list.addAll(temp);
	}
	
	private static boolean containsReadsAboveSize(ArrayList<Read> list, int size){
		for(Read r : list){
			if(r!=null && r.bases!=null){
				if(r.bases.length>size){
					assert(r.mate==null) : "Read of length "+r.bases.length+">"+size+". Paired input is incompatible with 'breaklength'";
					return true;
				}
			}
		}
		return false;
	}
	
	
	private long[] countReads(String fname1, String fname2, long maxReads){
		{
			String x=fname1.toLowerCase();
			if((x.equals("stdin") || x.startsWith("stdin.")) && !new File(fname1).exists()){
				throw new RuntimeException("Can't precount reads from standard in, only from a file.");
			}
		}
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, false, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Counting Reads");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		long count=0, count2=0, bases=0;
		
		while(reads!=null && reads.size()>0){
			count+=reads.size();
			for(Read r : reads){
				bases+=r.bases.length;
				count2++;
				if(r.mate!=null){
					bases+=r.mate.bases.length;
					count2++;
				}
			}
			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln, ln.list.isEmpty());
		errorState|=ReadWrite.closeStream(cris);
		return new long[] {count, count2, bases};
	}
	
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
		outstream.println("Other parameters and their defaults:\n");
		outstream.println("overwrite=false  \tOverwrites files that already exist");
		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
		outstream.println("fastawrap=80     \tLength of lines in fasta output");
		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	
	public void setSampleSeed(long seed){
		randy=new Random();
		if(seed>-1){randy.setSeed(seed);}
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;
	private String outsingle=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	/** Tracks names to ensure no duplicate names. */
	private final HashMap<String,Integer> nameMap1, nameMap2;
	private boolean uniqueNames=false;

	private boolean reverseComplimentMate=false;
	private boolean reverseCompliment=false;
	private boolean verifyinterleaving=false;
	private boolean trimBadSequence=false;
	private boolean deleteEmptyFiles=false;
	/** Add /1 and /2 to paired reads */
	private boolean addslash=false;

	private long maxReads=-1;
	private float samplerate=1f;
	private long sampleseed=-1;
	private boolean sampleReadsExact=false;
	private boolean sampleBasesExact=false;
	private long sampleReadsTarget=0;
	private long sampleBasesTarget=0;
	
	private boolean qtrimRight=false;
	private boolean qtrimLeft=false;
	private int forceTrimLeft=-1;
	private int forceTrimRight=-1;
	private byte trimq=4;
	private byte minAvgQuality=0;
	private int breakLength=0;
	private int maxReadLength=0;
	private int minReadLength=0;
	private float minLenFraction=0;
	/** Toss pair only if both reads are shorter than limit */ 
	private boolean requireBothBad=false;
	
	private boolean useSharedHeader;
	
	private boolean averageQualityByProbability=false;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	private final FileFormat ffoutsingle;
	
	private final boolean qtrim;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	private boolean colorspace=false;
	private boolean parsecustom=false;
	private boolean testsize=false;
	
	private Random randy;
	
}
