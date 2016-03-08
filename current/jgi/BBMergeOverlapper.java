package jgi;

import java.util.Arrays;

import kmer.KCountArray;
import stream.Read;
import align2.QualityTools;
import align2.Tools;
import dna.AminoAcid;

/**
 * @author Brian Bushnell
 * @date Apr 15, 2014
 *
 */
public final class BBMergeOverlapper {
	
	protected static final int mateByOverlap(Read a, Read b, int[] rvector, int minOverlap0, final int minOverlap, int margin,
			final int maxMismatches0, final int maxMismatches, final int minq) {
		minOverlap0=Tools.min(Tools.max(1, minOverlap0), minOverlap);
		assert(maxMismatches<=maxMismatches0);
		margin=Tools.max(margin, 0);
		assert(maxMismatches>=margin);
		if(rvector==null){rvector=new int[4];}
		if(USE_MAPPING){
			rvector[0]=100;
			rvector[1]=20;
			rvector[2]=0;
			//rvector[3]=20; Unused
			rvector[4]=0;
			return a.insertSizeMapped(IGNORE_MAPPING_STRAND);
		}
		
//		assert(false) : minOverlap0+", "+minOverlap+", "+margin+", "+maxMismatches0+", "+maxMismatches+", "+minq+", "+MIN_OVERLAP_INSERT;
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		
		int bestOverlap=-1;
//		int bestScore=-1; //This serves as the threshold for the minimum score to report. 
		int bestGood=-1;
		int bestBad=maxMismatches0;
//		float bestExpected=0;
		
		boolean ambig=false;
		final int maxOverlap=abases.length+bbases.length-Tools.max(minOverlap, MIN_OVERLAP_INSERT);
//		assert(false) : minOverlap+", "+maxOverlap;
//		System.err.print("\nm");
		
		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
//			System.err.print("\nn");
//			verbose=(insert==174);
			if(verbose){System.err.println("\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(abases.length+bbases.length-overlap));}
			
			
			int good=0, bad=0;
//			float expected=0;
			
			int istart=(overlap<=abases.length ? 0 : overlap-abases.length);
			int jstart=(overlap<=abases.length ? abases.length-overlap : 0);
//			System.err.print("o");
			
			for(int i=istart, j=jstart, badlim=bestBad+margin; i<overlap && i<bbases.length && j<abases.length && bad<badlim; i++, j++){
				assert(j>=0 && j<=abases.length && i>=0 && i<=bbases.length) : "\njstart="+jstart+", j="+j+
					", istart="+istart+", i="+i+" \n"+"overlap="+overlap+", a.length="+a.bases.length+
					", b.length="+b.bases.length+", bad="+bad+", badlim="+badlim+", good="+good;
				final byte ca=abases[j], cb=bbases[i];
				final byte qa=(aqual==null ? 20 : aqual[j]), qb=(bqual==null ? 20 : bqual[i]);
				byte q=QualityTools.qualsToPhred(qa, qb);
//				if(ca=='N' || cb=='N' || (aqual!=null && aqual[j]<minq) || (bqual!=null && bqual[i]<minq)){
				if(ca=='N' || cb=='N' || q<minq){
					//do nothing
				}else{
					assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
					if(ca==cb){good++;}
					else{bad++;}
//					expected+=QualityTools.qualsToProbError(qa, qb);
				}
			}
//			System.err.print("p");
			
//			System.err.println(overlap+", "+bestOverlap+", "+bestGood+", "+bestBad+", "+ambig);

//			System.err.print("a");
			if(bad*2<good){
				if(good>minOverlap){//Candidate
					if(bad<=bestBad){

						//					System.err.print("b");
						if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
//							ambig=(bestBad-bad<margin);
							if(bestBad-bad<margin){ambig=true;}
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
//							bestExpected=expected;
						}else if(bad==bestBad){
							ambig=true;
						}

						//					System.err.print("c");
						if(ambig && bestBad<margin){
							//						System.err.print("d");
							rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
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
					rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
					rvector[1]=bestGood;
					rvector[2]=bestBad;
					//				rvector[3]=bestSum;
					rvector[4]=(ambig ? 1 : 0);
					return -1;
				}
			}
//			System.err.print("h");
			
		}
//		System.err.println("i");
		
//		if(bestExpected<bestBad+.5f){ambig=true;}
//		if((bestExpected+.04f)*6f<bestBad){ambig=true;}
		
		if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}
		
		rvector[0]=((bestBad==0 ? GOOD_MULT_1 : GOOD_MULT_2)*bestGood-BAD_MULT*bestBad);
		rvector[1]=bestGood;
		rvector[2]=bestBad;
//		rvector[3]=bestSum;
		rvector[4]=(ambig ? 1 : 0);
//		assert(false) : Arrays.toString(rvector)+"\n"+new String(a.bases)+"\n"+new String(b.bases);
		
		return (bestOverlap<0 ? -1 : abases.length+bbases.length-bestOverlap);
	}
	
	
	/**
	 * TODO Use this
	 * @param a
	 * @param b
	 * @param overlap
	 * @return
	 */
	protected static final float expectedMismatches(Read a, Read b, int overlap) {
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		int good=0, bad=0;
		float expected=0;
		float goodWeighted=0;
		float badWeighted=0;
		
		if(aqual==null || bqual==null){return (overlap+0)/16;}

		int istart=(overlap<=abases.length ? 0 : overlap-abases.length);
		int jstart=(overlap<=abases.length ? abases.length-overlap : 0);
		//			System.err.print("o");

		for(int i=istart, j=jstart; i<overlap && i<bbases.length && j<abases.length; i++, j++){
			final byte ca=abases[j], cb=bbases[i];
			final byte qa=aqual[j], qb=bqual[i];
			
			if(ca=='N' || cb=='N'){
				//do nothing
			}else{
				assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) : (char)ca+", "+(char)cb;
				float prob=QualityTools.qualsToProbError(qa, qb);
				if(ca==cb){
					good++;
					goodWeighted+=(1-prob);
				}else{
					bad++;
					badWeighted+=prob;
				}
				expected+=prob;
			}
		}
		if(badWeighted>expected){return 0;}
		return expected;
	}
	
	private static final int countMismatches(Read a, Read b, int insert, int maxMismatches, byte minq){
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
				if(ca=='N' || cb=='N' || qa<minq || qb<minq){
					//do nothing
				}else{
					mismatches++;
					if(mismatches>maxMismatches){break;}
				}
			}
		}
		return mismatches;
	}
	
	protected static int minCoverage(final Read r, final KCountArray kca, final int k, boolean makeCanonical, int cutoff){
		final byte[] bases=r.bases;
		if(bases==null || bases.length<k){return cutoff;}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		
		int len=0;
		long kmer=0;
		int min=cutoff;
		
		for(int i=0; i<bases.length && min>=cutoff; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;

				if(len>=k){
					int cov=kca.read(kmer, k, makeCanonical);
					min=Tools.min(min, cov);
				}
			}
		}
		return min;
	}
	
	protected static int calcMinOverlapByEntropy(byte[] bases, int k, short[] counts, int minscore){
		final int bits=2*k;
		final int mask=~((-1)<<(bits));
		int kmer=0, len=0, ones=0, twos=0;
		
		Arrays.fill(counts, (short)0);
		
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			if(i<bases.length){
				final byte b=bases[j];
				if(!AminoAcid.isFullyDefined(b)){
					len=0;
					kmer=0;
				}else{
					len++;
					final int n=Dedupe.baseToNumber[b];
					kmer=((kmer<<2)|n)&mask;

					if(len>=k){
						counts[kmer]++;
						if(counts[kmer]==1){ones++;}
						else if(counts[kmer]==2){twos++;}
						if(ones*4+twos>=minscore){return i;}
					}
				}
			}
		}
		return bases.length+1;
	}
	
	public static int MIN_OVERLAP_INSERT=25;
	public static int BAD_MULT=6;
	public static int GOOD_MULT_1=8;
	public static int GOOD_MULT_2=400;

	/** Skip alignment and calculate insert from mapping info */ 
	public static boolean USE_MAPPING=false;
	public static final boolean IGNORE_MAPPING_STRAND=false;
	
	public static boolean verbose=false;
	
}
