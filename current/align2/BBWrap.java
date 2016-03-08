package align2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.Read;
import dna.Data;
import dna.Parser;
import dna.Timer;

/**
 * @author Brian Bushnell
 * @date Mar 27, 2014
 *
 */
public class BBWrap {
	
	public static void main(String[] args){
		BBWrap wrapper=new BBWrap();
		ArrayList<String> list=wrapper.parse(args);
		wrapper.execute(list);
	}

	private final ArrayList<String> parse(String[] args){

		sysout.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		sysout.println("BBMap version "+Shared.BBMAP_VERSION_STRING);
		
		if(Parser.parseHelp(args)){
//			printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer();
		t.start();
		
		Read.TO_UPPER_CASE=true;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
//			if("null".equalsIgnoreCase(b)){b=null;}
//			System.err.println("Processing "+arg);
			if(a.equals("path") || a.equals("root")){
				Data.setPath(b);
				args[i]=null;
			}else if(a.equals("mapper")){
				mapper=b;
				args[i]=null;
			}else if(a.equals("ref") || a.equals("reference") || a.equals("fasta")){
				ref=b;
				args[i]=null;
			}else if(a.equals("in") || a.equals("in1")){
				add(b, in1List);
				args[i]=null;
			}else if(a.equals("in2")){
				add(b, in2List);
				args[i]=null;
			}else if(a.equals("out") || a.equals("out1")){
				add(b, out1List);
				args[i]=null;
			}else if(a.equals("out2")){
				add(b, out2List);
				args[i]=null;
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("outmapped") || a.equals("outmapped1")){
				add(b, outm1List);
				args[i]=null;
			}else if(a.equals("outm2") || a.equals("outmapped2")){
				add(b, outm2List);
				args[i]=null;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("outunmapped") || a.equals("outunmapped1")){
				add(b, outu1List);
				args[i]=null;
			}else if(a.equals("outu2") || a.equals("outunmapped2")){
				add(b, outu2List);
				args[i]=null;
			}else if(a.equals("outb") || a.equals("outb1") || a.equals("outblack") || a.equals("outblack1") || a.equals("outblacklist") || a.equals("outblacklist1")){
				add(b, outb1List);
				args[i]=null;
			}else if(a.equals("outb2") || a.equals("outblack2") || a.equals("outblacklist2")){
				add(b, outb2List);
				args[i]=null;
			}else if(a.equals("qualityhistogram") || a.equals("qualityhist") || a.equals("qhist")){
				add(b, qhistList);
				args[i]=null;
			}else if(a.equals("matchhistogram") || a.equals("matchhist") || a.equals("mhist")){
				add(b, mhistList);
				args[i]=null;
			}else if(a.equals("inserthistogram") || a.equals("inserthist") || a.equals("ihist")){
				add(b, ihistList);
				args[i]=null;
			}else if(a.equals("bamscript") || a.equals("bs")){
				add(b, bsList);
				args[i]=null;
			}
		}
		
		ArrayList<String> list=new ArrayList<String>();
		for(String s : args){
			if(s!=null){
				list.add(s);
			}
		}
//		return list.toArray(new String[list.size()]);
		return list;
		
	}
	
	private void add(String s, ArrayList<String> list){
		if(s!=null && !"null".equals(s.toLowerCase())){
			String[] sa=s.split(",");
			for(String ss : sa){
				list.add(ss);
			}
		}
	}
	
	private void execute(ArrayList<String> base){
		for(int i=0; i<in1List.size(); i++){
			ArrayList<String> list=(ArrayList<String>) base.clone();
			
			if(i==0 && ref!=null){list.add("ref="+ref);}
			else if(i>0){list.add("indexloaded=t");}

			if(bsList.size()>0){list.add("bs="+bsList.get(i));}
			if(qhistList.size()>0){list.add("qhist="+qhistList.get(i));}
			if(mhistList.size()>0){list.add("mhist="+mhistList.get(i));}
			if(ihistList.size()>0){list.add("ihist="+ihistList.get(i));}
			if(in1List.size()>0){list.add("in="+in1List.get(i));}
			if(out1List.size()>0){list.add("out="+out1List.get(i));}
			if(outu1List.size()>0){list.add("outu="+outu1List.get(i));}
			if(outm1List.size()>0){list.add("outm="+outm1List.get(i));}
			if(outb1List.size()>0){list.add("outb="+outb1List.get(i));}
			if(in2List.size()>0){list.add("in2="+in2List.get(i));}
			if(out2List.size()>0){list.add("out2="+out2List.get(i));}
			if(outu2List.size()>0){list.add("outu2="+outu2List.get(i));}
			if(outm2List.size()>0){list.add("outm2="+outm2List.get(i));}
			if(outb2List.size()>0){list.add("outb2="+outb2List.get(i));}
			
			String[] args=list.toArray(new String[list.size()]);
			if(mapper==null || mapper.equalsIgnoreCase("bbmap")){
				BBMap.main(args);
			}else if(mapper.equalsIgnoreCase("bbmappacbio") || mapper.equalsIgnoreCase("pacbio")){
				BBMapPacBio.main(args);
			}else if(mapper.equalsIgnoreCase("bbmappacbioskimmer") || mapper.equalsIgnoreCase("pacbioskimmer") || mapper.equalsIgnoreCase("skimmer")){
				BBMapPacBioSkimmer.main(args);
			}else if(mapper.equalsIgnoreCase("bbmap5") || mapper.equalsIgnoreCase("5")){
				BBMap5.main(args);
			}else if(mapper.equalsIgnoreCase("bbmapacc") || mapper.equalsIgnoreCase("acc")){
				BBMapAcc.main(args);
			}else if(mapper.equalsIgnoreCase("bbsplit") || mapper.equalsIgnoreCase("bbsplitter")){
				BBSplitter.main(args);
			}
		}
	}

	private String ref;
	private String mapper="bbmap";

	private ArrayList<String> bsList=new ArrayList<String>();
	private ArrayList<String> qhistList=new ArrayList<String>();
	private ArrayList<String> mhistList=new ArrayList<String>();
	private ArrayList<String> ihistList=new ArrayList<String>();
	
	private ArrayList<String> in1List=new ArrayList<String>();
	private ArrayList<String> out1List=new ArrayList<String>();
	private ArrayList<String> outu1List=new ArrayList<String>();
	private ArrayList<String> outm1List=new ArrayList<String>();
	private ArrayList<String> outb1List=new ArrayList<String>();

	private ArrayList<String> in2List=new ArrayList<String>();
	private ArrayList<String> out2List=new ArrayList<String>();
	private ArrayList<String> outu2List=new ArrayList<String>();
	private ArrayList<String> outm2List=new ArrayList<String>();
	private ArrayList<String> outb2List=new ArrayList<String>();
	
	static PrintStream sysout=System.err;
	
}
