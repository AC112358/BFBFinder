package bfbf;

import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.Weights;
import main.java.bfbf.weights.FbrWeights;

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

    public boolean handle(PalindromeCollection collection, int l, double weight, double fbrWeight,
                          int prevFullSize, int prevSingletons) {
        boolean continueEnumeration;
        if (l == from - 1) {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, fbrWeight, w, fw, -1, minWeight,
                    minFbrWeight, sigCurves.get(0), sig, prevFullSize, prevSingletons);
            collection.unwrap();
        } else if (l < -1) {
            String s = collection.get(0).halfString();
            if (fw.getMiddleCharIndex() < 0 || s.charAt(s.length() - 1) == (char) ('A' + fw.getMiddleCharIndex())) {
                stream.println(s);
                ++stringCount;
            }
            continueEnumeration = stringCount < maxStrings;
        } else {
            collection.wrap();
            ;
            //System.out.println(collection);
            continueEnumeration = collection.allFoldings1(this, weight, fbrWeight, w, fw, l, minWeight,
                    minFbrWeight, sigCurves.get(l), sig, prevFullSize, prevSingletons);
            collection.unwrap();
        }
        return continueEnumeration;
    }

    public boolean handle(PalindromeCollection collection, int l, double weight) {
        boolean continueEnumeration;
        if (l == from - 1) {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, w, -1, minWeight, sigCurves.get(0), sig);
            collection.unwrap();
        } else if (l < -1) {
            String s = collection.get(0).halfString();
            stream.println(s);
            ++stringCount;
            continueEnumeration = stringCount < maxStrings;
        } else {
            collection.wrap();
            ;
            //System.out.println(collection);
            continueEnumeration = collection.allFoldings1(this, weight, w, l, minWeight, sigCurves.get(l), sig);
            collection.unwrap();
        }
        return continueEnumeration;
    }

    public int numOfStrings() {
        return stringCount;
    }


}