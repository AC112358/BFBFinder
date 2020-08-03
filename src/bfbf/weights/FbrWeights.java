package bfbf.weights;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FbrWeights extends Weights{

    protected FbrWeights(int length) {
        super(length);
    }

    public FbrWeights(int[] counts) {
        super(counts);
    }

    public FbrWeights(String inputStr) throws IllegalArgumentException {
        super(inputStr);
    }

    public FbrWeights(String inputStr, ErrorModel errorModel, double minWeight) throws IllegalArgumentException {
		super(inputStr, errorModel, minWeight);
    }

    public static void main(String[] args) {
        String str = "\n\n   3:0.5264,1.555  \n6:0.2,0.4,0.8";
        Weights w = new Weights(str);
    }

	public double getWeight(int ix, int fullSize, int prevFullSize, int currSingletons, int prevSingletons){
		int count = prevFullSize - fullSize + prevSingletons + currSingletons;
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
}
