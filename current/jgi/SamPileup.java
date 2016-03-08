package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

import stream.SamLine;

import align2.ReadStats;
import align2.Tools;

import dna.CoverageArray;
import dna.CoverageArray2;
import dna.CoverageArray3;
import dna.Data;
import dna.Gene;
import dna.Parser;
import dna.Scaffold;
import dna.Timer;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Jan 4, 2013
 *
 */
public class SamPileup {
	
	public static void main(String[] args){
		SamPileup sp=new SamPileup(args);

		Timer t=new Timer();
		t.start();
		
		sp.process();
		
		t.stop();
		Data.sysout.println();
		Data.sysout.println("Time: \t"+t);
		
	}
	
	public SamPileup(String[] args){
		for(String s : args){
			if(s.contains("=stdout")){Data.sysout=System.err;}
//			if(s.equals("in=stdin") || s.startsWith("in=stdin.")){SYSIN=true;}
		}
		System.err.println("Executing "+(this.getClass().getName())+" "+Arrays.toString(args)+"\n");
		
		boolean bs=false, setbs=false;
		boolean outset=false;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
//			System.err.println("Processing "+args[i]);
			
			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("ref") || a.equals("reference") || a.equals("fasta")){
				reference=b;
			}else if(a.equals("in") || a.equals("in1")){
				in=b;
			}else if(a.equals("out") || a.equals("outfile")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")){
//					System.err.println("No output file.");
					out=null;
				}else{
					out=b;
				}
				outset=true;
			}else if(a.equals("outsam") || a.equals("samout")){
				outsam=(b==null || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("outorf") || a.equals("orfout")){
				outorf=(b==null || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("orffasta") || a.equals("fastaorf")){
				orffasta=(b==null || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("basecov")){
				basecov=(b==null || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("bincov")){
				bincov=(b==null || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("delta")){
				deltaOnly=Tools.parseBoolean(b);
			}else if(a.equals("hist") || a.equals("histogram")){
				histogram=(b==null || b.equalsIgnoreCase("none")) ? null : b;
			}else if(a.equals("reads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("scafs") || a.equals("scaffolds")){
				initialScaffolds=Tools.max(128, (int)(Tools.min(Long.parseLong(b),2000000000)));
			}else if(a.equals("binsize")){
				binsize=Integer.parseInt(b);
			}else if(a.equals("32bit")){
				bits32=Tools.parseBoolean(b);
			}else if(a.equals("bitset") || a.equals("usebitset")){
				bs=Tools.parseBoolean(b);
				setbs=true;
			}else if(a.startsWith("nonzero") || a.equals("nzo")){
				NONZERO_ONLY=Tools.parseBoolean(b);
				System.err.println("Set NONZERO_ONLY to "+NONZERO_ONLY);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
				System.err.println("Set overwrite to "+overwrite);
			}else if(a.equalsIgnoreCase("twocolumn")){
				TWOCOLUMN=Tools.parseBoolean(b);
				System.err.println("Set TWOCOLUMN to "+TWOCOLUMN);
			}else if(a.equals("secondary") || a.equals("usesecondary")){
				USE_SECONDARY=Tools.parseBoolean(b);
				System.err.println("Set USE_SECONDARY_ALIGNMENTS to "+USE_SECONDARY);
			}else if(a.equals("keepshortbins") || a.equals("ksb")){
				KEEP_SHORT_BINS=Tools.parseBoolean(b);
				System.err.println("Set KEEP_SHORT_BINS to "+KEEP_SHORT_BINS);
			}else if(i>1){
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		if(setbs){
			USE_BITSETS=bs;
			USE_COVERAGE_ARRAYS=!bs;
			System.err.println("Set USE_BITSETS to "+USE_BITSETS);
		}else{
			if(histogram==null && basecov==null && bincov==null && outorf==null){//No need for coverage array!
				USE_COVERAGE_ARRAYS=false;
				if(TWOCOLUMN){//No need for bitset, either!
					USE_BITSETS=false;
				}else{
					USE_BITSETS=true;
				}
				System.err.println("Set USE_COVERAGE_ARRAYS to "+USE_COVERAGE_ARRAYS);
				System.err.println("Set USE_BITSETS to "+USE_BITSETS);
			}
		}
//		assert(false) : USE_COVERAGE_ARRAYS;
		
		if(maxReads<0){maxReads=Long.MAX_VALUE;}
		{
			final String a=(args.length>0 ? args[0] : null);
			final String b=(args.length>1 ? args[1] : null);
			if(in==null && a!=null && a.indexOf('=')<0 && (a.startsWith("stdin") || new File(a).exists())){in=a;}
			if(out==null && b!=null && b.indexOf('=')<0){out=b;}
			if(in==null){in="stdin";}
			if(out==null && !outset){
//				out="stdout";
//				System.err.println("Warning: output destination not set; producing no output.  To print to standard out, set 'out=stdout'");
				Data.sysout=System.err;
			}
		}
		assert(in!=null);
//		assert(out!=null || outset) : "Output file was not set.";
	}
	
	
	public void processOrfsFasta(String fname_in, String fname_out, HashMap<String, Scaffold> map){
		TextFile tf=new TextFile(fname_in, false, false);
		assert(!fname_in.equalsIgnoreCase(fname_out));
		TextStreamWriter tsw=new TextStreamWriter(fname_out, overwrite, false, true);
		tsw.start();

//		tsw.println("#refBases="+refBases);
		tsw.print("#mappedBases="+mappedBases+"\n");
		tsw.print("#name\tlength\tdepthSum\tavgDepth\tavgDepth/mappedBases\tminDepth\tmaxDepth\tmedianDepth\tstdDevDepth\tfractionCovered\n");
		
		String line;
		final StringBuilder sb=new StringBuilder(500);
		while((line=tf.nextLine())!=null){
			if(line.length()>1 && line.charAt(0)=='>'){
				
				String[] split=line.split(" # "); //' # ' used as delimiters
				
				String orfname=split[0].substring(1).trim(); //In case there are spaces around the ' # ' delimiters
				String scafname=orfname;
				if(scafname.contains("_")){//PRODIGAL pads _1 to the name of the first orf of a scaffold, and etc
					int last=scafname.lastIndexOf('_');
					boolean numeric=false;
					for(int i=last+1; i<scafname.length(); i++){
						if(Character.isDigit(scafname.charAt(i))){numeric=true;}
						else{numeric=false; break;}
					}
					if(numeric){scafname=scafname.substring(0, last);}
				}
				
				int start=Integer.parseInt(split[1].trim());
				int stop=Integer.parseInt(split[2].trim());
				int strand=Integer.parseInt(split[3].trim());
				if(strand==1){strand=Gene.PLUS;}else{strand=Gene.MINUS;}
				Orf orf=new Orf(orfname, start, stop, (byte)strand);
				
				Scaffold scaf=map.get(scafname);
//				if(scaf==null){scaf=map.get(orfname);}
				
//				assert(scaf!=null) : "\nCan't find scaffold for ("+orf+")\nfrom line\n"+line+"\n";
//				assert(orf.start>=0 && orf.stop<scaf.length) : "\norf goes out of scaffold bounds.\n"+orf+"\n"+scaf+"\n";
				
				if(scaf==null){
					System.err.println("Can't find scaffold for ("+orf+")\nfrom line\n"+line+"\nscafname='"+scafname+"'\norfname='"+orfname+"'");
					if(ABORT_ON_ERROR){
						tsw.poison();
						throw new RuntimeException("Aborting.");
					}
				}
				if(orf.start<0 && orf.stop>=scaf.length){
					Data.sysout.println("orf goes out of scaffold bounds.\n"+orf+"\n"+scaf);
					if(ABORT_ON_ERROR){
						tsw.poison();
						throw new RuntimeException("Aborting.");
					}
				}

				if(scaf!=null){
					CoverageArray ca=(CoverageArray)scaf.obj;
					orf.readCoverageArray(ca);
				}

				sb.append(orf.name).append('\t');
				sb.append(orf.length()).append('\t');
				sb.append(orf.baseDepth).append('\t');
				sb.append(String.format("%.4f", orf.avgCoverage())).append('\t');
				sb.append(orf.avgCoverage()/mappedBases);

				sb.append('\t');
				sb.append(orf.minDepth).append('\t');
				sb.append(orf.maxDepth).append('\t');
				sb.append(orf.medianDepth).append('\t');
				sb.append(String.format("%.4f",orf.stdevDepth)).append('\t');
				sb.append(String.format("%.4f",orf.fractionCovered()));

				sb.append('\n');
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		
		tsw.poison();
		tsw.waitForFinish();
	}
	
	
	public void process(){
		refBases=0;
		mappedBases=0;
		ArrayList<Scaffold> list=new ArrayList<Scaffold>(initialScaffolds);
		HashMap<String, Scaffold> table=new HashMap<String, Scaffold>(initialScaffolds);
		TextFile tf=new TextFile(in, false, false);
		String line=null;
		
		String program=null;
		String version=null;
		
		boolean bbmap=false;
		float bbversion=-1;
		
		final TextStreamWriter tsw=(outsam==null ? null : new TextStreamWriter(outsam, overwrite, false, true));
		if(outsam!=null){tsw.start();}
		
		for(line=tf.nextLine(); line!=null && line.startsWith("@"); line=tf.nextLine()){
			if(tsw!=null){tsw.println(line);}
			
			final String[] split=line.split("\t");
			final String a=split[0];
			
			if(a.equals("@SQ")){
				Scaffold sc=new Scaffold(split);
				if(COUNT_GC){sc.basecount=new long[6];}
				assert(!table.containsKey(sc.name)) : "\nDuplicate scaffold name!\n"+sc+"\n\n"+table.get(sc.name);
				table.put(sc.name, sc);
				list.add(sc);
				refBases+=sc.length;
//				sc.obj=new CoverageArray2(table.size(), sc.length+1);
//				Data.sysout.println("Made scaffold "+sc.name+" of length "+sc.length);
			}else if(a.equals("@PG")){
				for(String s : split){
					if(s.startsWith("PN:")){
						String s2=s.substring(3);
						if(s2.equalsIgnoreCase("bbmap") || s2.startsWith("BBMap")){bbmap=true;}
						if(program==null){program=Data.forceIntern(s.substring(3));}
					}else if(s.startsWith("VN:")){
						if(bbmap && bbversion<0){bbversion=Float.parseFloat(s.substring(3));}
						if(version==null){version=Data.forceIntern(s.substring(3));}
					}
				}
			}else if(a.equals("@RG")){
				//Do nothing
			}else if(a.equals("@HD")){
				//Do nothing
			}else if(a.equals("@CO")){
				//Do nothing
			}else{
//				assert(false) : line;
			}
		}
		
		if(reference!=null){
			TextFile tf2=new TextFile(reference, false, false);
			Scaffold sc=null;
			int len=0;
			final long[] acgtn=new long[6];
			for(String s=tf2.nextLine(); s!=null; s=tf2.nextLine()){
				if(s.startsWith(">")){
					if(sc!=null){
						sc.length=len;
						sc.gc=(float)((acgtn[1]+acgtn[2])*1d/Data.max(1, acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3]));
						sc=null;
						len=0;
						Arrays.fill(acgtn, 0);
					}
					
					String name=s.substring(1);
					sc=table.get(name);
					if(ADD_FROM_REF && sc==null){
						sc=new Scaffold(name, 0);
						System.err.println("Warning - SAM header did not include "+name);
						table.put(name, sc);
					}
				}else{
					len+=s.length();
					for(int i=0; i<s.length(); i++){
						acgtn[charToNum[s.charAt(i)]]++;
					}
				}
			}
			if(sc!=null){
				sc.length=len;
				sc.gc=(float)((acgtn[1]+acgtn[2])*1d/Data.max(1, acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3]));
				sc=null;
				len=0;
				Arrays.fill(acgtn, 0);
			}
		}
		
		boolean err=false;
		for(; line!=null; line=tf.nextLine()){
			if(tsw!=null){tsw.println(line);}
			
			if(line.length()==0){
				
			}else if(line.charAt(0)=='@'){
				if(!err){
					System.err.println("Unexpected header line: "+line);
					System.err.println("This should not cause problems, and is probably due to concatenated sam files.\n" +
							"Supressing future unexpected header warnings.");
					err=true;
				}
				
				if(line.startsWith("@SQ")){
					String[] split=line.split("\t");
					Scaffold sc=new Scaffold(split);
//					if(bbmap && bbversion<=17 && sc.length>1000){sc.length-=1000;}
					if(!table.containsKey(sc.name)){
						if(COUNT_GC){sc.basecount=new long[6];}
						table.put(sc.name, sc);
						list.add(sc);
						refBases+=sc.length;
					}
				}
			}else{

				SamLine sl=new SamLine(line);
				if(sl.mapped() && (USE_SECONDARY || sl.primary())){
					mappedBases+=sl.seq.length;
					final Scaffold scaf=table.get(new String(sl.rname()));
					assert(scaf!=null) : "Can't find "+new String(sl.rname());
					final int a=Tools.max(sl.start(), 0);
					final int b=Tools.min(sl.stop2(), scaf.length-1);
					scaf.basehits+=(b-a+1);
					
					if(USE_COVERAGE_ARRAYS){
						if(scaf.obj==null){
							scaf.obj=(bits32 ? new CoverageArray3(table.size(), scaf.length+1) : new CoverageArray2(table.size(), scaf.length+1));
						}
						CoverageArray ca=(CoverageArray)scaf.obj;
						ca.incrementRange(a, b, 1);
					}else if(USE_BITSETS){
						if(scaf.obj==null){
							scaf.obj=new BitSet(scaf.length+1);
						}
						BitSet bs=(BitSet)scaf.obj;
						bs.set(a, b+1);
					}
//					assert(false) : a+", "+b+", "+scaf.length;
					if(sl.seq!=null && scaf.basecount!=null){
						final long[] counts=scaf.basecount;
//						final String seq=sl.seq;
//						for(int i=0; i<seq.length(); i++){
//							counts[charToNum[seq.charAt(i)]]++;
//						}
						final byte[] seq=sl.seq;
						for(int i=0; i<seq.length; i++){
							counts[charToNum[seq[i]]]++;
						}
					}
				}
			}
		}
		tf.close();
		if(tsw!=null){tsw.poison();}
		
//		OutputStream os=ReadWrite.getOutputStream(out, false);
//		PrintWriter pw=new PrintWriter(os);
		
//		for()
		
		final TextStreamWriter tsw2=(out==null ? null : new TextStreamWriter(out, overwrite, false, true));
		
		if(tsw2!=null){
			tsw2.start();
			if(TWOCOLUMN){
				tsw2.println("ID\tAvg_fold");
			}else if(COUNT_GC){
				tsw2.println("ID\tAvg_fold\tLength\tRef_GC\tBase_Coverage\tRead_GC");
//				tsw2.println("ID\tLength\tRef_GC\tAvg_fold\tBase_Coverage\tRead_GC");
			}else{
				tsw2.println("ID\tAvg_fold\tLength\tRef_GC\tBase_Coverage");
//				tsw2.println("ID\tLength\tRef_GC\tAvg_fold\tBase_Coverage");
			}
		}
		
		totalScaffolds=list.size();
		if(USE_COVERAGE_ARRAYS || USE_BITSETS /*histogram!=null || tsw2!=null*/){
			final long[] hist=(USE_COVERAGE_ARRAYS ? new long[Character.MAX_VALUE+1] : null);

			for(Scaffold scaf : list){
				final long sum=scaf.basehits;
				int covered=0;
				if(USE_COVERAGE_ARRAYS){
					CoverageArray ca=(CoverageArray)scaf.obj;
					if(ca!=null){
						for(int i=0; i<scaf.length; i++){
							int x=ca.get(i);
							hist[x]++;
//							sum+=x;
							if(x>0){covered++;}
						}
					}
				}else if(USE_BITSETS){
//					sum+=scaf.basehits;
					BitSet bs=(BitSet)scaf.obj;
					covered=(bs==null ? 0 : bs.cardinality());
				}
				
				if(sum>0){
					scaffoldsWithCoverage++;
				}
				//			pw.print(scaf.name);
				if(tsw2!=null && (sum>0 || !NONZERO_ONLY)){
					if(TWOCOLUMN){
						tsw2.print(String.format("%s\t%.4f\n", scaf.name, sum/(double)scaf.length));
					}else if(COUNT_GC){
						long[] bc=scaf.basecount;
						double gc=(bc[1]+bc[2])*1d/Data.max(1, bc[0]+bc[1]+bc[2]+bc[3]);
						tsw2.print(String.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%.4f\n", scaf.name, sum/(double)scaf.length, scaf.length, scaf.gc, covered*100d/scaf.length, gc));
					}else{
						tsw2.print(String.format("%s\t%.4f\t%d\t%.4f\t%.4f\n", scaf.name, sum/(double)scaf.length, scaf.length, scaf.gc, covered*100d/scaf.length));
					}
				}
				totalCoveredBases+=covered;
			}
			
			if(tsw2!=null){tsw2.poison();}

			if(histogram!=null && hist!=null){writeHist(histogram, hist);}
		}
		
		if(basecov!=null){writeCoveragePerBase(basecov, list, deltaOnly);}
		if(bincov!=null){
			if(KEEP_SHORT_BINS){
				writeCoveragePerBaseBinned2(bincov, list, binsize);
			}else{
				writeCoveragePerBaseBinned(bincov, list, binsize);
			}
		}
		
		if(orffasta!=null){
			processOrfsFasta(orffasta, outorf, table);
		}
		
		if(tsw2!=null){tsw2.waitForFinish();}
		if(tsw!=null){tsw.waitForFinish();}
		

		double depthCovered=mappedBases*1.0/refBases;
		double pctScaffoldsWithCoverage=scaffoldsWithCoverage*100.0/totalScaffolds;
		double pctCovered=totalCoveredBases*100.0/refBases;
		
		Data.sysout.println(String.format("\nAverage coverage:                    \t%.2f", depthCovered));
		Data.sysout.println(String.format("Percent scaffolds with any coverage: \t%.2f", pctScaffoldsWithCoverage));
		if(USE_COVERAGE_ARRAYS || USE_BITSETS){
			Data.sysout.println(String.format("Percent of reference bases covered:  \t%.2f", pctCovered));
		}
	}
	
	public static void writeHist(String fname, long[] counts){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#Coverage\tnumBases\n");
		int max=0;
		for(max=counts.length-1; max>0 && counts[max]==0; max--){}
		for(int i=0; i<=max; i++){
			long x=counts[i];
			tsw.print(i+"\t"+x+"\n");
		}
		tsw.poison();
		tsw.waitForFinish();
	}
	
	public static void writeCoveragePerBase(String fname, ArrayList<Scaffold> list, boolean deltaOnly){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, true);
		tsw.start();
		tsw.print("#RefName\tPos\tCoverage\n");
		
		for(Scaffold scaf : list){
			int last=-1;
			CoverageArray ca=(CoverageArray)scaf.obj;
			for(int i=0; i<scaf.length; i++){
				int x=(ca==null ? 0 : ca.get(i));
				if(!deltaOnly || x!=last){
					tsw.print(scaf.name+"\t"+(i+1)+"\t"+x+"\n");
					last=x;
				}
			}
		}
		
		tsw.poison();
//		tsw.waitForFinish();
	}
	
	/** Note.  As written, this will truncate all trailing bases of each scaffold's length modulo binsize.
	 * For example, with binsize 1000, the last 500 bases of a 1500 base scaffold will be ignored. 
	 * @param fname
	 * @param list
	 * @param binsize
	 */
	public static void writeCoveragePerBaseBinned(String fname, ArrayList<Scaffold> list, int binsize){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#RefName\tCov\tPos\tRunningPos\n");
		
		long running=0;
		final float invbin=1f/binsize;
		for(Scaffold scaf : list){
			if(scaf.length>=binsize){
				CoverageArray ca=(CoverageArray)scaf.obj;
				int lastPos=-1, nextPos=binsize-1;
				long sum=0;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					sum+=x;
					if(i>=nextPos){
//						float bin=(i-lastPos);
//						tsw.print(String.format("%s\t%.1f\t%d\t%d\n", scaf.name, sum/bin, (i+1), running));
						tsw.print(String.format("%s\t%.2f\t%d\t%d\n", scaf.name, sum*invbin, (i+1), running));
						lastPos=i;
						running+=binsize;
						nextPos+=binsize;
						sum=0;
					}
				}
			}
		}
		
		tsw.poison();
		tsw.waitForFinish();
	}
	
	/** This version will NOT truncate all trailing bases of each scaffold's length modulo binsize.
	 * @param fname
	 * @param list
	 * @param binsize
	 */
	public static void writeCoveragePerBaseBinned2(String fname, ArrayList<Scaffold> list, int binsize){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#RefName\tCov\tPos\tRunningPos\n");
		
		long running=0;
		for(Scaffold scaf : list){
			CoverageArray ca=(CoverageArray)scaf.obj;
			int lastPos=-1, nextPos=binsize-1;
			long sum=0;
			final int lim=scaf.length-1;
			for(int i=0; i<scaf.length; i++){
				int x=(ca==null ? 0 : ca.get(i));
				sum+=x;
				if(i>=nextPos || i==lim){
					int bin=(i-lastPos);
					tsw.print(String.format("%s\t%.1f\t%d\t%d\n", scaf.name, sum/(float)bin, (i+1), running));
					running+=bin;
					nextPos+=binsize;
					lastPos=i;
					sum=0;
				}
			}
		}
		
		tsw.poison();
		tsw.waitForFinish();
	}
	
	
	public static boolean COUNT_GC=true;
	public static boolean verbose=false; 
	/** Permission to overwrite existing files */
	public static boolean overwrite=false;
	/** Permission to append to existing files */
	public static boolean append=false;
	public static boolean TWOCOLUMN=false;
	public static boolean ADD_FROM_REF=false;
	public static boolean NONZERO_ONLY=false;
	public static boolean USE_COVERAGE_ARRAYS=true;
	public static boolean USE_BITSETS=false;
	public static boolean deltaOnly=false;
	/** Process secondary alignments */
	public static boolean USE_SECONDARY=true;
	public static boolean ABORT_ON_ERROR=true;
	public static boolean KEEP_SHORT_BINS=true;

	public long maxReads=-1;
	public int initialScaffolds=4096;
	public String in=null;
	public String out=null;
	public String outsam=null;
	public String outorf=null;
	public String reference=null;
	public String histogram=null;
	public String basecov=null;
	public String bincov=null;
	public String orffasta=null;
	public boolean bits32=false;
	public int binsize=1000;

	public long refBases=0;
	public long mappedBases=0;
	public long totalCoveredBases=0;
	public long scaffoldsWithCoverage=0;
	public long totalScaffolds=0;
	
	private static final byte[] charToNum=makeCharToNum();
	private static byte[] makeCharToNum() {
		byte[] r=new byte[256];
		Arrays.fill(r, (byte)4);
		r['a']=r['A']=0;
		r['c']=r['C']=1;
		r['g']=r['G']=2;
		r['t']=r['T']=3;
		r['\n']=r['\r']=r['>']=r['@']=r['+']=5;
		return r;
	}
	
	
	
	
}
