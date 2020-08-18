package bfbf;


import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.FbrWeights;
import bfbf.weights.Weights;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class AllPairwiseBFBStringAccumulator extends FoldingHandler {

    protected int maxStrings;
    protected int stringCount;
    protected List<Set<String>[]> pairwiseStrings;
    protected List<StringBuilder[]> stringBuilders;

    public AllPairwiseBFBStringAccumulator(Weights w, int from, int to, double minWeight) {
        this(w, from, to, minWeight, -1);
    }

    public AllPairwiseBFBStringAccumulator(Weights w, FbrWeights fw, int from, int to,
                                           double minWeight, double minFbrWeight) {
        super(w, fw, from, to, minWeight, 0);
    }

    public AllPairwiseBFBStringAccumulator(Weights w, int from, int to, double minWeight,
                                           int maxStrings) {
        super(w, from, to, minWeight);
        if (maxStrings < 0) {
            this.maxStrings = Integer.MAX_VALUE;
        } else {
            this.maxStrings = maxStrings;
        }

        stringCount = 0;

        pairwiseStrings = new ArrayList<>();
        stringBuilders = new ArrayList<>();

        pairwiseStrings.add(null);
        stringBuilders.add(null);

        ensureLength(to);
    }

    public void ensureLength(int k) {
        for (int i = pairwiseStrings.size(); i < k; ++i) {
            Set<String>[] sets = new Set[i];
            pairwiseStrings.add(sets);
            StringBuilder[] builders = new StringBuilder[i];
            stringBuilders.add(builders);
            for (int j = 0; j < i; ++j) {
                sets[j] = new HashSet<>();
                builders[j] = new StringBuilder();
            }
        }
    }

    public void clear() {
        stringCount = 0;
        for (int i = 1; i < pairwiseStrings.size(); ++i) {
            Set<String>[] sets = pairwiseStrings.get(i);
            //			StringBuilder[] builders = stringBuilders.get(i);
            for (int j = 0; j < i; ++j) {
                sets[j].clear();
                //				builders[j].setLength(0);
            }
        }
    }

    public boolean handle(PalindromeCollection collection, int l, double weight1, double weight2,
                          int prevFullSize, int prevSingletons) {
        return handle(collection, l, weight1);
    }

    public boolean handle(PalindromeCollection collection, int l, double weight) {
        boolean continueEnumeration;
        if (l == from - 1) {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, w, l, minWeight, sigCurves.get(0), sig);
            collection.unwrap();
        } else if (l < -1) {
            int[] seq = collection.get(0).seq(from);
            int length = seq.length / 2;

            for (int i = from + 1; i < to; ++i) {
                for (int j = from; j < i; ++j) {
                    stringBuilders.get(i)[j].setLength(0);
                }
            }

            for (int q = 0; q < length; ++q) {
                char currChar = (char) ('A' + seq[q]);
                for (int p = from; p < seq[q]; ++p) {
                    stringBuilders.get(seq[q])[p].append(currChar);
                }
                for (int p = seq[q] + 1; p < to; ++p) {
                    stringBuilders.get(p)[seq[q]].append(currChar);
                }
            }

            for (int i = from + 1; i < to; ++i) {
                for (int j = from; j < i; ++j) {
                    pairwiseStrings.get(i)[j].add(stringBuilders.get(i)[j].toString());
                }
            }

            ++stringCount;
            continueEnumeration = stringCount < maxStrings;
        } else {
            collection.wrap();
            continueEnumeration = collection.allFoldings1(this, weight, w, l, minWeight, sigCurves.get(l - from), sig);
            collection.unwrap();
        }

        return continueEnumeration;
    }

    public int getStringCount() {
        return stringCount;
    }

    public void printAll() {
        printAll(System.out);
    }

    public void printAll(PrintStream out) {
        char c1, c2;
        for (int i = from + 1; i < to; ++i) {
            Set<String>[] sets = pairwiseStrings.get(i);
            c2 = (char) ('A' + i);
            for (int j = from; j < i; ++j) {
                if (!sets[j].isEmpty()) {
                    c1 = (char) ('A' + j);
                    out.println("Projections over " + c1 + ", " + c2 + " (" +
                            sets[j].size() + "):");
                    for (String str : sets[j]) {
                        out.println(str);
                    }
                    out.println();
                }
            }
        }
    }

    public void setFromTo(int from, int to) {
        this.from = from;
        this.to = to;
    }

}
