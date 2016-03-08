package jgi;

import stream.Read;
import align2.Tools;
import dna.AminoAcid;

/**
 * @author Brian Bushnell
 * @date Apr 15, 2014
 *
 */
public final class BBMergeOverlapper {
	
	public static final int mateByOverlap(Read a, Read b, int[] rvector, final int minOverlap0, final int minOverlap) {
		if(rvector==null){rvector=new int[4];}
		if(USE_MAPPING){
			rvector[0]=100;
			rvector[1]=20;
			rvector[2]=0;
			//rvector[3]=20; Unused
			rvector[4]=0;
			return a.insertSizeMapped(IGNORE_MAPPING_STRAND);
		}
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		
		int bestOverlap=-1;
//		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=DEFAULT_BADLIMIT;
		final int margin=2;
		
		boolean ambig=false;
		final int maxOverlap=abases.length+bbases.length-Tools.max(minOverlap, MIN_OVERLAP_INSERT);
//		assert(false) : minOverlap+", "+maxOverlap;
//		System.err.print("\nm");
		
		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
//			System.err.print("\nn");
//			verbose=(insert==174);
			if(verbose){System.err.println("\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(abases.length+bbases.length-overlap));}
			
			
			int tested=0;
			int good=0, bad=0;
			
			int istart=(overlap<=abases.length ? 0 : overlap-abases.length);
			int jstart=(overlap<=abases.length ? abases.length-overlap : 0);
//			System.err.print("o");
			
			for(int i=istart, j=jstart, badlim=bestBad+margin; i<overlap && i<bbases.length && j<abases.length && bad<badlim; i++, j++){
				assert(j>=0 && j<=abases.length && i>=0 && i<=bbases.length) : "\njstart="+jstart+", j="+j+", istart="+istart+", i="+i+" \n"+
						"overlap="+overlap+", a.length="+a.bases.length+", b.length="+b.bases.length+", bad="+bad+", badlim="+badlim+", good="+good+", tested="+tested;
				byte ca=abases[j], cb=bbases[i];
				if(ca=='N' || cb=='N' || (aqual!=null && aqual[j]<MIN_QUALITY) || (bqual!=null && bqual[i]<MIN_QUALITY)){
					//do nothing
				}else{
					assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
					tested++;
					if(ca==cb){good++;}
					else{bad++;}
				}
			}
//			System.err.print("p");
			
//			System.err.println(overlap+", "+bestOverlap+", "+bestGood+", "+bestBad+", "+ambig);

//			System.err.print("a");
			if(good>minOverlap){//Candidate
				if(bad<=bestBad){

//					System.err.print("b");
					if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
						if(bad>bestBad-margin){ambig=true;}
						bestOverlap=overlap;
						bestBad=bad;
						bestGood=good;
					}else if(bad==bestBad){
						ambig=true;
					}

//					System.err.print("c");
					if(ambig && bestBad<margin){
//						System.err.print("d");
						rvector[0]=((bestBad==0 ? 8 : 4)*bestGood-6*bestBad);
						rvector[1]=bestGood;
						rvector[2]=bestBad;
//						rvector[3]=bestSum;
						rvector[4]=(ambig ? 1 : 0);
//						System.err.print("e");
						return -1;
					}

//					System.err.print("f");
				}
			}else if(bad<margin){
//				System.err.print("g");
				ambig=true;
				rvector[0]=((bestBad==0 ? 8 : 4)*bestGood-6*bestBad);
				rvector[1]=bestGood;
				rvector[2]=bestBad;
//				rvector[3]=bestSum;
				rvector[4]=(ambig ? 1 : 0);
				return -1;
			}
//			System.err.print("h");
			
		}
//		System.err.println("i");
		
		rvector[0]=((bestBad==0 ? 8 : 4)*bestGood-6*bestBad);
		rvector[1]=bestGood;
		rvector[2]=bestBad;
//		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
		
		return (bestOverlap<0 ? -1 : abases.length+bbases.length-bestOverlap);
	}
	
	
	private static final int countMismatches(Read a, Read b, int insert, int maxMismatches){
		final int lengthSum=a.bases.length+b.bases.length;
		if(insert>=lengthSum){return 0;}
		final int overlap=Tools.min(insert, lengthSum-insert);
		
		int mismatches=0;
		
		
		int start1=(insert>a.bases.length ? a.bases.length-overlap : 0);
		int start2=(insert>=b.bases.length ? 0 : b.bases.length-overlap);
		
		while(start1<0 || start2<0){start1++; start2++;}
		for(int i=start1, j=start2; i<a.bases.length && j<b.bases.length; i++, j++){
			final byte ca=a.bases[i], cb=b.bases[j];
			if(ca!=cb){
				final byte qa=a.quality[i], qb=b.quality[j];
				if(ca=='N' || cb=='N' || qa<7 || qb<7){
					//do nothing
				}else{
					mismatches++;
					if(mismatches>maxMismatches){break;}
				}
			}
		}
		return mismatches;
	}
	
	public static int DEFAULT_BADLIMIT=3;
	public static int MIN_OVERLAP_INSERT=16;
	public static byte MIN_QUALITY=7;

	/** Skip alignment and calculate insert from mapping info */ 
	public static boolean USE_MAPPING=false;
	public static final boolean IGNORE_MAPPING_STRAND=false;
	
	public static boolean verbose=false;
	
}
