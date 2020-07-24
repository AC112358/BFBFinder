package bfbf;

import java.util.List;


/**
 * An abstract class for objects implementing a generic "handle" function
 * over palindrome collections.  
 * 
 * 
 * @author zakov
 *
 */
public abstract class FoldingHandler {
	
	protected Weights w;
	protected int from, to; // segment interval begin and end indices
	protected List<SigCurve> sigCruves;
	protected Signature sig;
	protected double minWeight;
	

	public FoldingHandler(Weights w, int from, int to, double minWeight) {
		super();
		sig = new Signature();
		set(w, from, to, minWeight);
	}


	public void set(Weights w, int from, int to, double minWeight) {
		this.w = w;
		this.from = from;
		this.to = to;
		sigCruves = SigCurve.sigCurves(w, minWeight, from, to);
		this.minWeight = minWeight;
	}


	/**
	 * An abstract procedure applied over a BFB palindrome collection.
	 * 
	 * @param collection the input collection.
	 * @param l base layer level.
	 * @param weight the collection's weight.
	 * 
	 * @return true if the procedure succeeded, false otherwise.
	 */
	public abstract boolean handle(PalindromeCollection collection, int l, double weight);
}
