package dna;

import java.io.File;

import jgi.CalcTrueQuality;

import stream.ConcurrentGenericReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import fileIO.ByteFile;
import fileIO.ReadWrite;

/**
 * @author Brian Bushnell
 * @date Mar 21, 2014
 *
 */
public class Parser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public Parser(){}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean parse(String arg, String a, String b){
		if(isJavaFlag(arg)){return true;}
		
		if(parseZip(arg, a, b)){return true;}
		if(parseSam(arg, a, b)){return true;}
		if(parseFasta(arg, a, b)){return true;}
		if(parseCommonStatic(arg, a, b)){return true;}

		if(parseFiles(arg, a, b)){return true;}
		if(parseCommon(arg, a, b)){return true;}
		if(parseQuality(arg, a, b)){return true;}
		if(parseTrim(arg, a, b)){return true;}
		if(parseInterleaved(arg, a, b)){return true;}
		
		return false;
	}

	public boolean parseCommon(String arg, String a, String b){
		if(a.equals("reads") || a.equals("maxreads")){
			maxReads=Long.parseLong(b);
		}else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
			assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
		}else if(a.equals("sampleseed")){
			sampleseed=Long.parseLong(b);
		}else if(a.equals("t") || a.equals("threads")){
			Shared.THREADS=Tools.max(Integer.parseInt(b), 1);
		}else if(a.equals("overwrite") || a.equals("ow")){
			overwrite=Tools.parseBoolean(b);
		}else if(a.equals("testsize")){
			testsize=Tools.parseBoolean(b);
		}else if(a.equals("breaklen") || a.equals("breaklength") || a.equals("maxlength") || a.equals("maxlen")){
			breakLength=Integer.parseInt(b);
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseQuality(String arg, String a, String b){
		if(a.equals("ignorebadquality") || a.equals("ibq")){
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
		}else if(a.equals("fakequality") || a.equals("qfake")){
			FASTQ.FAKE_QUAL=Byte.parseByte(b);
		}else if(a.equals("qauto")){
			FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseInterleaved(String arg, String a, String b){
		if(a.equals("testinterleaved")){
			FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
			System.err.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
			setInterleaved=true;
		}else if(a.equals("forceinterleaved")){
			FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
			System.err.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			setInterleaved=true;
		}else if(a.equals("interleaved") || a.equals("int")){
			if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
			else{
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				setInterleaved=true;
			}
		}else if(a.equals("overrideinterleaved")){
			boolean x=Tools.parseBoolean(b);
			ReadStreamByteWriter.ignorePairAssertions=x;
			if(x){setInterleaved=true;}
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseTrim(String arg, String a, String b){
		if(a.equals("ftl") || a.equals("forcetrimleft")){
			forceTrimLeft=Integer.parseInt(b);
		}else if(a.equals("ftr") || a.equals("forcetrimright")){
			forceTrimRight=Integer.parseInt(b);
		}else if(a.equals("trim") || a.equals("qtrim")){
			if(b==null){qtrimRight=qtrimLeft=true;}
			else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrimLeft=true;qtrimRight=false;}
			else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrimLeft=false;qtrimRight=true;}
			else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrimLeft=qtrimRight=true;}
			else{qtrimRight=qtrimLeft=Tools.parseBoolean(b);}
		}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
			if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
				TrimRead.optimalMode=true;
				TrimRead.optimalBias=Float.parseFloat(b);
				assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
			}else{
				TrimRead.optimalMode=Tools.parseBoolean(b);
			}
		}else if(a.equals("trimright") || a.equals("qtrimright")){
			qtrimRight=Tools.parseBoolean(b);
		}else if(a.equals("trimleft") || a.equals("qtrimleft")){
			qtrimLeft=Tools.parseBoolean(b);
		}else if(a.equals("trimq") || a.equals("trimquality")){
			trimq=Byte.parseByte(b);
		}else if(a.equals("requirebothbad") || a.equals("rbb")){
			requireBothBad=Tools.parseBoolean(b);
		}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
			minReadLength=Integer.parseInt(b);
		}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
			minLenFraction=Float.parseFloat(b);
		}else if(a.equals("minavgquality") || a.equals("maq")){
			minAvgQuality=Byte.parseByte(b);
		}else if(a.equals("trimBadSequence") || a.equals("tbs")){
			trimBadSequence=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}
	
	public boolean parseFiles(String arg, String a, String b){
		if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
			in1=b;
		}else if(a.equals("in2") || a.equals("input2")){
			in2=b;
		}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
			out1=b;
		}else if(a.equals("out2") || a.equals("output2")){
			out2=b;
		}else if(a.equals("qfin") || a.equals("qfin1")){
			qfin1=b;
		}else if(a.equals("qfout") || a.equals("qfout1")){
			qfout1=b;
		}else if(a.equals("qfin2")){
			qfin2=b;
		}else if(a.equals("qfout2")){
			qfout2=b;
		}else if(a.equals("extin")){
			extin=b;
		}else if(a.equals("extout")){
			extout=b;
		}else{
			return false;
		}
		return true;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean parseCommonStatic(String arg, String a, String b){
		if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription")){
			Shared.TRIM_READ_COMMENTS=Tools.parseBoolean(b);
		}else if(a.equals("tuc") || a.equals("touppercase")){
			Read.TO_UPPER_CASE=Tools.parseBoolean(b);
		}else if(a.equals("lctn") || a.equals("lowercaseton")){
			Read.LOWER_CASE_TO_N=Tools.parseBoolean(b);
		}else if(a.equals("tossbrokenreads") || a.equals("tbr")){
			boolean x=Tools.parseBoolean(b);
			Read.NULLIFY_BROKEN_QUALITY=x;
			ConcurrentGenericReadInputStream.REMOVE_DISCARDED_READS=x;
		}else if(a.equals("build") || a.equals("genome")){
			Data.setGenome(Integer.parseInt(b));
		}else if(a.equals("bf1")){
			ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
			ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
		}else if(a.equals("bf2")){
			ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
			ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
		}else{
			return false;
		}
		return true;
	}

	public static boolean parseZip(String arg, String a, String b){
		if(a.equals("ziplevel") || a.equals("zl")){
			int x=Integer.parseInt(b);
			if(x>=0){
				ReadWrite.ZIPLEVEL=Tools.min(x, 9);
			}
		}else if(a.equals("usegzip") || a.equals("gzip")){
			ReadWrite.USE_GZIP=Tools.parseBoolean(b);
		}else if(a.equals("usepigz") || a.equals("pigz")){
			if(b!=null && Character.isDigit(b.charAt(0))){
				int zt=Integer.parseInt(b);
				if(zt<1){ReadWrite.USE_PIGZ=false;}
				else{
					ReadWrite.USE_PIGZ=true;
					if(zt>1){
						ReadWrite.MAX_ZIP_THREADS=zt;
						ReadWrite.ZIP_THREAD_DIVISOR=1;
					}
				}
			}else{ReadWrite.USE_PIGZ=Tools.parseBoolean(b);}
		}else if(a.equals("usegunzip") || a.equals("gunzip")){
			ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
		}else if(a.equals("useunpigz") || a.equals("unpigz")){
			ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}

	public static boolean parseSam(String arg, String a, String b){
		if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
			SamLine.VERSION=Float.parseFloat(b);
		}else if(a.equals("mdtag") || a.equals("md")){
			SamLine.MAKE_MD_TAG=Tools.parseBoolean(b);
		}else if(a.equals("xstag") || a.equals("xs")){
			SamLine.MAKE_XS_TAG=true;
			if(b!=null){
				b=b.toLowerCase();
				if(b.startsWith("fr-")){b=b.substring(3);}
				if(b.equals("ss") || b.equals("secondstrand")){
					SamLine.XS_SECONDSTRAND=true;
				}else if(b.equals("fs") || b.equals("firststrand")){
					SamLine.XS_SECONDSTRAND=false;
				}else if(b.equals("us") || b.equals("unstranded")){
					SamLine.XS_SECONDSTRAND=false;
				}else{
					SamLine.MAKE_XS_TAG=Tools.parseBoolean(b);
				}
			}
			SamLine.setxs=true;
		}else if(a.equals("intronlen") || a.equals("intronlength")){
			SamLine.INTRON_LIMIT=Integer.parseInt(b);
			SamLine.setintron=true;
		}else if(a.equals("idtag")){
			SamLine.MAKE_IDENTITY_TAG=Tools.parseBoolean(b);
		}else if(a.equals("xmtag") || a.equals("xm")){
			SamLine.MAKE_XM_TAG=Tools.parseBoolean(b);
		}else if(a.equals("stoptag")){
			SamLine.MAKE_STOP_TAG=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}

	public static boolean parseFasta(String arg, String a, String b){
		if(a.equals("fastareadlen") || a.equals("fastareadlength")){
			FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
			FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
		}else if(a.equals("fastaminread") || a.equals("fastaminlen") || a.equals("fastaminlength")){
			FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
		}else if(a.equals("fastawrap")){
			FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
		}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
			int x=Integer.parseInt(b);
			FastaReadInputStream.MIN_READ_LEN=(x>0 ? x : Integer.MAX_VALUE);
		}else{
			return false;
		}
		return true;
	}
	
	public static boolean parseQualityAdjust(String arg, String a, String b){
		if(a.equals("q102matrix") || a.equals("q102m")){
			CalcTrueQuality.q102matrix=b;
		}else if(a.equals("qbpmatrix") || a.equals("bqpm")){
			CalcTrueQuality.qbpmatrix=b;
		}else if(a.equals("loadq102")){
			CalcTrueQuality.q102=Tools.parseBoolean(b);
		}else if(a.equals("loadqbp")){
			CalcTrueQuality.qbp=Tools.parseBoolean(b);
		}else if(a.equals("loadq10")){
			CalcTrueQuality.q10=Tools.parseBoolean(b);
		}else if(a.equals("loadq12")){
			CalcTrueQuality.q12=Tools.parseBoolean(b);
		}else if(a.equals("loadqb012")){
			CalcTrueQuality.qb012=Tools.parseBoolean(b);
		}else if(a.equals("loadqb234")){
			CalcTrueQuality.qb234=Tools.parseBoolean(b);
		}else if(a.equals("loadqp")){
			CalcTrueQuality.qp=Tools.parseBoolean(b);
		}else if(a.equals("adjustquality") || a.equals("adjq")){
			TrimRead.ADJUST_QUALITY=Tools.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}

	public static boolean isJavaFlag(String arg){
		if(arg==null){return false;}
		if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.startsWith("-Xmn") || arg.equals("-ea") || arg.equals("-da")){return true;}
		if(arg.startsWith("Xmx") || arg.startsWith("Xms") || arg.startsWith("Xmn")){
			return arg.length()>3 && Character.isDigit(arg.charAt(3));
		}
		return false;
	}
	

	/** Return true if the user seems confused */
	public static boolean parseHelp(String[] args){
		if(args==null || args.length==0 || (args.length==1 && args[0]==null)){return true;}
		if(args.length>1){return false;}
		final String s=args[0].toLowerCase();
		return s.equals("-h") || s.equals("-help") || s.equals("--help") 
				|| s.equals("-version") || s.equals("--version") || s.equals("?") || s.equals("-?") || (s.equals("help") && !new File(s).exists());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int forceTrimLeft=-1;
	public int forceTrimRight=-1;

	public long maxReads=-1;
	public float samplerate=1f;
	public long sampleseed=-1;

	public boolean qtrimLeft=false;
	public boolean qtrimRight=false;
	
	public byte trimq=6;
	public byte minAvgQuality=0;
	public int minReadLength=0;
	public float minLenFraction=0;
	public int breakLength=0;
	/** Toss pair only if both reads are shorter than limit */ 
	public boolean requireBothBad=false;
	public boolean trimBadSequence=false;

	public boolean overwrite=false;
	public boolean testsize=false;
	
	public boolean setInterleaved=false;
	
	public String in1=null;
	public String in2=null;
	
	public String qfin1=null;
	public String qfin2=null;

	public String out1=null;
	public String out2=null;

	public String qfout1=null;
	public String qfout2=null;
	
	public String extin=null;
	public String extout=null;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static byte qin=-1;
	public static byte qout=-1;
	
}
