package align2;

/**
 * @author Brian Bushnell
 * @date Aug 5, 2013
 *
 */
public abstract class BandedAligner {
	
	public BandedAligner(int width_){
		maxWidth=Tools.max(width_, 3)|1;
		assert(maxWidth>=3) : "width<3 : "+width_+" -> "+maxWidth;
		assert(big>maxWidth/2);
	}
	
	public static final BandedAligner makeBandedAligner(int width_){
		BandedAligner ba=(Shared.USE_JNI ? new BandedAlignerJNI(width_) : new BandedAlignerConcrete(width_));
		return ba;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignForward(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignForwardRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignReverse(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignReverseRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	protected void fillBig(int[] array){
		final int lim=array.length-1;
		for(int i=1; i<lim; i++){array[i]=big;}
	}
	
	/** Score is lastRow-edits */
	public final int score(){
		return lastRow-lastEdits+1;
	}
	
	/** Position of min value in array (meaning the best alignment) relative to the middle of the array. */
	protected int lastOffset(int[] array, int halfWidth){
		final int center=halfWidth+1;
		int minLoc=center;
		for(int i=1; i<=halfWidth; i++){
			if(array[center+i]<array[minLoc]){minLoc=center+i;}
			if(array[center-i]<array[minLoc]){minLoc=center-i;}
		}
		return center-minLoc;
	}
	
	protected int penalizeOffCenter(int[] array, int halfWidth){
		final int center=halfWidth+1;
		int edits=array[center];
		for(int i=1; i<=halfWidth; i++){
			array[center+i]=Tools.min(big, array[center+i]+i);
			edits=Tools.min(edits, array[center+i]);
			array[center-i]=Tools.min(big, array[center-i]+i);
			edits=Tools.min(edits, array[center-i]);
		}
		return edits;
	}
	
	/** Final row aligned in last alignment. */
	public int lastRow;
	/** Final edits value in last alignment. */
	public int lastEdits;

	/** Position of min value in array (meaning the best alignment) relative to the middle of the array.
	 * Positive value is to the right (ref sequence longer than query), negative value left (ref shorter than query) */
	protected int lastOffset;
	
	public int lastRefLoc;
	public int lastQueryLoc;
	
	public final int maxWidth;
	
	protected static final int big=999;
	public static boolean verbose=false;
	/** Penalizes non-length-neutral alignments.  
	 * This Causes query-to-ref alignment to yield same score as ref-to-query alignment, which is useful for assertions.  */ 
	public static boolean penalizeOffCenter=true;
	
}
