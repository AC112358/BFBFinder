package bfbf;

import bfbf.weights.Weights;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SigCurve {
    protected List<Signature> sigs;
    protected TDoubleList weights;
    private double minWeight;

    protected SigCurve(double minWeight) {
        sigs = new ArrayList<>(1); // a series of lexicographically increasing signatures
        weights = new TDoubleArrayList(1); // a series of decreasing weights
        sigs.add(new Signature(1));
        weights.add(1);
        this.minWeight = minWeight;
    }

    protected SigCurve(SigCurve prev, Weights w, int l, Signature s) {
        sigs = new ArrayList<>(); // a series of lexicographically increasing signatures
        weights = new TDoubleArrayList(); // a series of decreasing weights
        this.minWeight = prev.minWeight;

        int maxCount = w.getMaxCount(l);
        for (int currCount = w.getMinCount(l); currCount <= maxCount; ++currCount) {
            double currWeight = w.getWeight(l, currCount);
            int maxIx = prev.findSmallerWeightIx(minWeight / currWeight);
            int minIx = prev.findMaxSignature(currCount, maxIx);
            for (int j = minIx; j < maxIx; ++j) {
                //System.out.println(prev.sigs);
                s.setTo(prev.sigs.get(j));
                if (s.minDecrement(currCount)) {
                    double weight = currWeight * prev.weights.get(j);
                    int insertionIx = Collections.binarySearch(sigs, s);
                    if (insertionIx < 0) {
                        insertionIx = -insertionIx - 1;
                        if (insertionIx < sigs.size() && weights.get(insertionIx) > weight) {
                            // the next signature has a higher weight
                            continue;
                        } else {
                            sigs.add(insertionIx, new Signature(s));
                            weights.insert(insertionIx, weight);
                        }
                    } else {
                        if (weight > weights.get(insertionIx)) {
                            weights.set(insertionIx, weight);
                        }
                    }
                    for (int p = insertionIx - 1; p >= 0 && weights.get(p) <= weight; --p) {
                        sigs.remove(p);
                        weights.removeAt(p);
                    }
                }
            }
        }
    }

    public static List<SigCurve> sigCurves(Weights weights, double minWeight) {
        return sigCurves(weights, minWeight, 0, weights.length());
    }

    public static List<SigCurve> sigCurves(Weights weights, double minWeight, int start) {
        return sigCurves(weights, minWeight, start, weights.length());
    }


    public static List<SigCurve> sigCurves(Weights weights, double minWeight, int start, int end) {
        int k = end - start;
        Signature s = new Signature();
        List<SigCurve> curves = new ArrayList<>(k);
        SigCurve currCurve, prevCurve = new SigCurve(minWeight);
        curves.add(prevCurve);

        for (int l = start; l < end; ++l) {
            currCurve = new SigCurve(prevCurve, weights, l, s);
            curves.add(currCurve);
            prevCurve = currCurve;
        }

        return curves;
    }

    private int findMaxSignature(int count, int max) {
        int min = 0;
        while (min < max) {
            int mid = (max + min) / 2;
            if (sigs.get(mid).canUnfold(count)) {
                max = mid;
            } else {
                min = mid + 1;
            }
        }
        return max;
    }

    /**
     * @param weight
     * @return the first index such that all weights as from this index are
     * strictly smaller than the input weight
     */
    private int findSmallerWeightIx(double weight) {
        int max = weights.size();
        int min = 0;
        while (min < max) {
            int mid = (max + min) / 2;
            if (weights.get(mid) < weight) {
                max = mid;
            } else {
                min = mid + 1;
            }
        }
        return max;
    }

    public boolean isEmpty() {
        return sigs.isEmpty();
    }

    public boolean withinValidRange(Signature s, double weight) {
        int rank = Collections.binarySearch(sigs, s);
        if (rank < 0) {
            rank = -rank - 1;
        }

        return (rank < sigs.size() &&
                weight * weights.get(rank) >= minWeight);
    }

    public boolean withinValidRange(Signature sig) {
        return sig.compareTo(sigs.get(sigs.size() - 1)) <= 0;
    }

    public int size() {
        return sigs.size();
    }

    public void add(int insertionIx, Signature signature, double weight) {
        sigs.add(insertionIx, signature);
        weights.insert(insertionIx, weight);
    }
}
