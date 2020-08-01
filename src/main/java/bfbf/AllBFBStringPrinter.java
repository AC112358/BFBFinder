package bfbf;

import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.Weights;

import java.io.PrintStream;

public class AllBFBStringPrinter extends FoldingHandler {

    protected PrintStream stream;
    protected int maxStrings;
    int stringCount;

    public AllBFBStringPrinter(Weights w, int from, int to, double minWeight,
                               PrintStream stream) {
        this(w, from, to, minWeight, stream, -1);
    }

    public AllBFBStringPrinter(Weights w, int from, int to, double minWeight,
                               PrintStream stream, int maxStrings) {
        super(w, from, to, minWeight);
        this.stream = stream;
        if (maxStrings < 0) {
            this.maxStrings = Integer.MAX_VALUE;
        } else {
            this.maxStrings = maxStrings;
        }
        stringCount = 0;
    }

    public boolean handle(PalindromeCollection collection, int l, double weight) {
        boolean continueEnumeration;
        if (l == from - 1) {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, w, -1, minWeight, sigCurves.get(0), sig);
            collection.unwrap();
        } else if (l < -1) {
            stream.println(collection.get(0).halfString());
            ++stringCount;
            continueEnumeration = stringCount < maxStrings;
        } else {
            collection.wrap();
            ;
            continueEnumeration = collection.allFoldings1(this, weight, w, l, minWeight, sigCurves.get(l), sig);
            collection.unwrap();
        }
        return continueEnumeration;
    }

    public int numOfStrings() {
        return stringCount;
    }


}
