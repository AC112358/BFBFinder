package bfbf;

import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.Weights;

public class AllBFBStringCounter extends FoldingHandler {
	
	private int[] countBox;
	private int maxCount;

	public AllBFBStringCounter(Weights w, int from, int to, double minWeight, int[] countBox) {
		this(w, from, to, minWeight, countBox, -1);
	}

	
	public AllBFBStringCounter(Weights w, int from, int to, double minWeight, int[] countBox, int maxCount) {
		super(w, from, to, minWeight);
		if (maxCount < 0){
			this.maxCount = Integer.MAX_VALUE;
		}
		else{
			this.maxCount = maxCount;
		}
		this.countBox = countBox;
	}

	public boolean handle(PalindromeCollection collection, int l, double weight) {
		boolean continueEnumeration;
		
		if (l == from-1){
			collection.wrap();
			continueEnumeration = collection.allFoldings1(this, weight, w, -1, minWeight, sigCurves.get(0), sig);
			collection.unwrap();
		}
		else if (l < -1){
			++countBox[0];
			continueEnumeration = countBox[0] < maxCount; 
			if (countBox[0] % 1000 == 0){
				System.out.println(countBox[0]/1000 + "K");
			}
		}
		else{
			collection.wrap(); //sig.accommodate(collection)
			continueEnumeration = collection.allFoldings1(this, weight, w, l, minWeight, sigCurves.get(l), sig);
			collection.unwrap();
		}
		return continueEnumeration;
	}


}
