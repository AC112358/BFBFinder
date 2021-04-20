package bfbf.weights;

import java.util.regex.Pattern;

public class FbrPosNegWeights extends FbrWeights {
	//protected String input;
	//protected ErrorModel error;
	protected int numLoopsIndex;

	protected static final Pattern countVecPtrn = Pattern.compile(
			"\\[\\s*((\\d+|-1)(\\s*,\\s*(\\d+|-1))*)?\\s*\\]");

	protected FbrWeights wPos;
	protected FbrWeights wNeg;


    public FbrPosNegWeights(String inputStr, ErrorModel errorModel, double minWeight) throws IllegalArgumentException {
		super(1);
    	int endIndex = inputStr.indexOf("][");
		String inputStrPos = inputStr.substring(0, endIndex + 1);
		String inputStrNeg = inputStr.substring(endIndex + 1);
    	wPos = new FbrWeights(inputStrPos, errorModel, minWeight, countVecPtrn);
    	wNeg = new FbrWeights(inputStrNeg, errorModel, minWeight, countVecPtrn);
		numLoopsIndex = -1;
		int[] posCounts = wPos.getCounts();
		int[] negCounts = wNeg.getCounts();
		int[] totCounts = new int[posCounts.length];
		for (int i = 0; i < totCounts.length; i++){
			totCounts[i] = posCounts[i] + negCounts[i];
		}
		super.updateCounts(totCounts, errorModel, minWeight);
    	//error = errorModel;
		//input = inputStr;
    }

    public void setNumLoopsIndex(int newIndex){
    	numLoopsIndex = newIndex;
	}

	public void updateCounts(int[] countsPos, int[] countsNeg, int incrementIndex,
							 ErrorModel errorModel,
							 double minWeight){
    	wPos.updateCounts(countsPos, incrementIndex, errorModel, minWeight);
		wNeg.updateCounts(countsNeg, incrementIndex, errorModel, minWeight);
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
		int fbrCount = prevFullSize - fullSize + prevSingletons + currSingletons;
    	//System.out.println("fullSize = " + fullSize + ", prevFullSize = " + prevFullSize +
		//		", prevSingletons = " + prevSingletons + ", currSingletons = " + currSingletons);
    	int fbrCountNeg = prevFullSize - (fullSize - currSingletons);
    	int fbrCountPos = fbrCount - fbrCountNeg;
    	//System.out.println("fbrCountNeg = " + fbrCountNeg);
    	//System.out.println("fbrCountPos = " + fbrCountPos);
    	//System.out.println("fbrCount = " + fbrCount);
		//System.out.println("foldback count for " + ix + " equals " + count);
		return wPos.getWeight(ix, fbrCountPos) * wNeg.getWeight(ix, fbrCountNeg);
	}
/*
	public int getMinCount(int l, double w, int fullSize, int prevFullSize, int prevSingletons){
    	int minFbrCountPos = wPos.getMinCount(l, w);
		int minFbrCountNeg = wNeg.getMinCount(l, w);
		int count = Math.min(minFbrCountNeg, minFbrCountPos);
		int minSingletons = count - prevFullSize + fullSize - prevSingletons;
		return minSingletons;
	}

	public int getMaxCount(int l, double w, int fullSize, int prevFullSize, int prevSingletons) {
		int maxFbrCountPos = wPos.getMaxCount(l, w);
		int maxFbrCountNeg = wNeg.getMaxCount(l, w);
		int count = Math.max(maxFbrCountNeg, maxFbrCountPos);
		int maxSingletons = count - prevFullSize + fullSize - prevSingletons;
		return maxSingletons;
	}
*/
	public int[] getPosCounts(){
    	return wPos.getCounts();
	}
	public int[] getNegCounts(){
		return wNeg.getCounts();
	}



}
