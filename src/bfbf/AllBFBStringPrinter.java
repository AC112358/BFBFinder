package bfbf;

import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.Weights;
import bfbf.weights.FbrWeights;

import java.io.PrintStream;

public class AllBFBStringPrinter extends FoldingHandler {

    protected PrintStream stream;
    protected int maxStrings;
    int stringCount;

    public AllBFBStringPrinter(Weights w, FbrWeights fw, int from, int to, double minWeight, double minFbrWeight,
                               PrintStream stream) {
        this(w, fw, from, to, minWeight, minFbrWeight, stream, -1);
    }

    public AllBFBStringPrinter(Weights w, FbrWeights fw, int from, int to, double minWeight, double minFbrWeights,
                               PrintStream stream, int maxStrings) {
        super(w, fw, from, to, minWeight, minFbrWeights);
        this.stream = stream;
        if (maxStrings < 0) {
            this.maxStrings = Integer.MAX_VALUE;
        } else {
            this.maxStrings = maxStrings;
        }
        stringCount = 0;
    }

    public boolean handle(PalindromeCollection collection, int l, double weight, double fbrWeight) {
        System.out.println("handle with two args has been called");
        boolean continueEnumeration;
        if (l == from - 1) {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, fbrWeight, w, fw, -1, minWeight, minFbrWeight, sigCurves.get(0), sig);
            collection.unwrap();
        } else if (l < -1) {
            stream.println(collection.get(0).halfString());
            ++stringCount;
            continueEnumeration = stringCount < maxStrings;
        } else {
            collection.wrap();
            ;
            System.out.println(collection);
            continueEnumeration = collection.allFoldings1(this, weight, fbrWeight, w, fw, l, minWeight, minFbrWeight, sigCurves.get(l), sig);
            collection.unwrap();
        }
        return continueEnumeration;
    }

    public boolean handle(PalindromeCollection collection, int l, double weight) {
        boolean continueEnumeration;
        if (l == from - 1) {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, fbrWeight, w, -1, minWeight, sigCurves.get(0), sig);
            collection.unwrap();
        } else if (l < -1) {
            stream.println(collection.get(0).halfString());
            ++stringCount;
            continueEnumeration = stringCount < maxStrings;
        } else {
            collection.wrap();
            ;
            System.out.println(collection);
            continueEnumeration = collection.allFoldings1(this, weight, w, l, minWeight, sigCurves.get(l), sig);
            collection.unwrap();
        }
        return continueEnumeration;
    }

    public int numOfStrings() {
        return stringCount;
    }


}
