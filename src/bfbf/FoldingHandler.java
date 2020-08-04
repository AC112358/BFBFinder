package bfbf;

import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.Weights;
import bfbf.weights.FbrWeights;

import java.util.List;


/**
 * An abstract class for objects implementing a generic "handle" function
 * over palindrome collections.
 *
 * @author zakov
 */
public abstract class FoldingHandler {

    protected Weights w;
    protected FbrWeights fw;
    protected int from, to; // segment interval begin and end indices
    protected List<SigCurve> sigCurves;
    protected Signature sig;
    protected double minWeight;
    protected double minFbrWeight;


    public FoldingHandler(Weights w, FbrWeights fw, int from, int to, double minWeight, double minFbrWeight) {
        super();
        sig = new Signature();
        set(w, fw, from, to, minWeight, minFbrWeight);
    }

    public FoldingHandler(Weights w, int from, int to, double minWeight) {
        super();
        sig = new Signature();
        set(w, null, from, to, minWeight, 0);
    }


    public void set(Weights w, FbrWeights fw, int from, int to, double minWeight, double minFbrWeight) {
        this.w = w;
        this.fw = fw;
        this.from = from;
        this.to = to;
        sigCurves = SigCurve.sigCurves(w, minWeight, from, to);
        this.minWeight = minWeight;
        this.minFbrWeight = minFbrWeight;
    }


    /**
     * An abstract procedure applied over a BFB palindrome collection.
     *
     * @param collection the input collection.
     * @param l          base layer level.
     * @param weight1     the collection's count weight.
     * @param weight2     the collection's fbr weight
     * @return true if the procedure succeeded, false otherwise.
     */
    public abstract boolean handle(PalindromeCollection collection, int l, double weight1, double weight2);

    /**
     * An abstract procedure applied over a BFB palindrome collection.
     *
     * @param collection the input collection.
     * @param l          base layer level.
     * @param weight     the collection's count weight.
     * @return true if the procedure succeeded, false otherwise.
     */
    public abstract boolean handle(PalindromeCollection collection, int l, double weight);
}
