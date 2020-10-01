package bfbf.weights;

import java.util.regex.Pattern;

public class FbrWeights extends Weights {
	//protected String input;
	//protected ErrorModel error;
	protected int numLoopsIndex;

	protected static final Pattern countVecPtrn = Pattern.compile(
			"\\[\\s*((\\d+|-1)(\\s*,\\s*(\\d+|-1))*)?\\s*\\]");

    protected FbrWeights(int length) {
        super(length);
        numLoopsIndex = -1;
    }

    public FbrWeights(int[] counts) {
        super(counts);
		numLoopsIndex = -1;
    }

    public FbrWeights(String inputStr) throws IllegalArgumentException {
        super(inputStr);
		numLoopsIndex = -1;
    }

    public FbrWeights(String inputStr, ErrorModel errorModel, double minWeight) throws IllegalArgumentException {
		super(inputStr, errorModel, minWeight, countVecPtrn);
		numLoopsIndex = -1;
    	//error = errorModel;
		//input = inputStr;
    }

    public void setNumLoopsIndex(int newIndex){
    	numLoopsIndex = newIndex;
	}

	public void updateCounts(int[] counts, int incrementIndex, ErrorModel errorModel, double minWeight){
    	Weights w = errorModel.getWeights(counts, minWeight);
    	this.minCounts = w.minCounts;
    	this.weights = w.weights;
    	this.heaviestCounts = w.heaviestCounts;
		numLoopsIndex = incrementIndex;
	}
/*
    public String getInput(){
    	return input;
	}

	public void setInput(String s){
    	input = s;
	}

	public ErrorModel getError(){
    	return error;
	}
*/
    public static void main(String[] args) {
        String str = "\n\n   3:0.5264,1.555  \n6:0.2,0.4,0.8";
        Weights w = new Weights(str);
    }

    public int getNumLoopsIndex(){
    	return numLoopsIndex;
	}

	public double getWeight(int ix, int fullSize, int prevFullSize, int currSingletons, int prevSingletons){
    	if (prevFullSize < 0){
    		return 1;
		}
		int count = prevFullSize - fullSize + prevSingletons + currSingletons;
		//System.out.println("foldback count for " + ix + " equals " + count);
		return getWeight(ix, count);
	}

	public int getMinCount(int l, double w, int fullSize, int prevFullSize, int prevSingletons){
    	int count = getMinCount(l, w);
		int minSingletons = count - prevFullSize + fullSize - prevSingletons;
		return minSingletons;
	}

	public int getMaxCount(int l, double w, int fullSize, int prevFullSize, int prevSingletons){
    	int count = getMaxCount(l, w);
		int maxSingletons = count - prevFullSize + fullSize - prevSingletons;
		return maxSingletons;
	}

	public int getMinCount(int l, double w) {
		if (l < 0) return 1;
		int c = heaviestCounts[l];
		for (; c >= 1 && getWeight(l, c - 1) >= w; --c) ;
		return c;
	}

	public double getWeight(int ix, int count) {
    	if (heaviestCounts[ix] < 0){
    		return 1;
		}
    	return super.getWeight(ix, 2 * (count/2));
	}


}
