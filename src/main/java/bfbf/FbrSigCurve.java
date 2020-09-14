package bfbf;

import bfbf.weights.Weights;
import bfbf.weights.FbrWeights;
import bfbf.distances.Distances;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class FbrSigCurve{
    protected List<nmSignaturePair> nmSigs;
    //protected List<Signature> sigs;
    protected TDoubleList weights;
    protected TDoubleList fbrWeights;
    //protected List<Integer> nmSums;
    protected List<Integer> counts;
    private double minWeight;
    private double minFbrWeight;
    protected List<Integer> nmStreaks;
    protected List<Integer> nmRegions;

    protected FbrSigCurve(double minWeight, double minFbrWeight) {
        nmSigs = new ArrayList<>(1);
        nmSignaturePair firstElmt = new nmSignaturePair(1, new Signature(1));
        nmSigs.add(firstElmt);
        weights = new TDoubleArrayList(1); // a series of decreasing weights
        fbrWeights = new TDoubleArrayList(1); //a series of decreasing fbr weights
        counts = new ArrayList<>(1);
        nmStreaks = new ArrayList<>(1);
        nmRegions = new ArrayList<>(1);
        weights.add(1);
        fbrWeights.add(1);
        counts.add(1);
        nmStreaks.add(1);
        nmRegions.add(1);
        this.minWeight = minWeight;
        this.minFbrWeight = minFbrWeight;
    }


    private int findSmallerWeightIx(double weight, TDoubleList weightsList) {
        int max = weightsList.size();
        int min = 0;
        while (min < max) {
            int mid = (max + min) / 2;
            if (weightsList.get(mid) < weight) {
                max = mid;
            } else {
                min = mid + 1;
            }
        }
        return max;
    }

    private int findSmallerWeightIx(double weight, int min, int max) {
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

    private int findMaxSignature(int count, int min, int max) {
        while (min < max) {
            int mid = (max + min) / 2;
            if (nmSigs.get(mid).sig.canUnfold(count)) {
                max = mid;
            } else {
                min = mid + 1;
            }
        }
        return max;
    }

    private int findMaxSignature(int count, int max) {
        int min = 0;
        while (min < max) {
            int mid = (max + min) / 2;
            if (nmSigs.get(mid).sig.canUnfold(count)) {
                max = mid;
            } else {
                min = mid + 1;
            }
        }
        return max;
    }

    private FbrSigCurve(FbrSigCurve prev, Weights w, FbrWeights fw, int l, Signature s, Distances dist) {
        nmSigs = new ArrayList<>(); // a series of lexicographically increasing signatures
        weights = new TDoubleArrayList(); // a series of decreasing weights
        fbrWeights = new TDoubleArrayList();
        counts = new ArrayList<>();
        this.minWeight = prev.minWeight;
        this.minFbrWeight = prev.minFbrWeight;
        nmStreaks = new ArrayList<>();
        nmRegions = new ArrayList<>();

        int maxCount = w.getMaxCount(l);
        for (int currCount = w.getMinCount(l); currCount <= maxCount; ++currCount) {
            int numElmts = 0;
            double currWeight = w.getWeight(l, currCount);
            for (int i = 0; i < prev.nmStreaks.size(); i++) {
                int maxIx = prev.findSmallerWeightIx(minWeight / currWeight, numElmts, numElmts + prev.nmStreaks.get(i));
                int minIx = prev.findMaxSignature(currCount, numElmts, maxIx);
                for (int j = minIx; j < maxIx; ++j) {
                    s.setTo(prev.nmSigs.get(j).sig);
                    int prevNmSum = prev.nmSigs.get(j).nmSum;
                    int prevN = 0;
                    if (prevNmSum > 0) {
                        prevN = prev.counts.get(j);
                    }
                    int prevM = prevNmSum - prevN;
                    double weight = currWeight * prev.weights.get(j);
                    int minFbrCount = fw.getMinCount(l, minFbrWeight, prevN, currCount, prevM);
                    int maxFbrCount = fw.getMaxCount(l, minFbrWeight, prevN, currCount, prevM);
                    for (int currEpsilon = minFbrCount; currEpsilon <= maxFbrCount; ++currEpsilon) {
                        s.setTo(prev.nmSigs.get(j).sig);
                        if (!s.minTwoDecrements(currCount, currCount - currEpsilon)) {
                            continue;
                        }
                        double currFbrWeight = fw.getWeight(l, prevN, prevM, currCount, currEpsilon);
                        double fbrWeight = currFbrWeight * prev.fbrWeights.get(j);
                        nmSignaturePair searchKey = new nmSignaturePair(currCount + currEpsilon, s);
                        //int insertionIx = Collections.binarySearch(sigs, s);
                        int nmIx = Collections.binarySearch(nmRegions,currCount + currEpsilon);
                        int insertionIx = Collections.binarySearch(nmSigs, searchKey);
                        int tempNmSum = currCount + currEpsilon;
                        if (nmIx < 0) {
                            insertionIx = -insertionIx - 1;
                            nmSigs.add(insertionIx, new nmSignaturePair(tempNmSum,
                                    new Signature(s)));
                            weights.insert(insertionIx, weight);
                            fbrWeights.insert(insertionIx, fbrWeight);
                            counts.add(insertionIx, currCount);
                            nmIx = -nmIx - 1;
                            if (nmIx >= nmStreaks.size()){
                                nmStreaks.add(1);
                                nmRegions.add(tempNmSum);
                            } else {
                                nmStreaks.add(nmIx, 1);
                                nmRegions.add(nmIx, tempNmSum);
                            }
                        } else if (insertionIx < 0) {
                            insertionIx = -insertionIx - 1;
                            if (insertionIx < nmSigs.size() &&
                                    dist.firstPairDominates(weights.get(insertionIx), fbrWeights.get(insertionIx),
                                            weight, fbrWeight)) {
                                // the next signature has a higher weight
                                continue;
                            } else {
                                nmSigs.add(insertionIx, new nmSignaturePair(currCount + currEpsilon,
                                        new Signature(s)));
                                weights.insert(insertionIx, weight);
                                fbrWeights.insert(insertionIx, fbrWeight);
                                counts.add(insertionIx, currCount);
                                nmStreaks.set(nmIx, nmStreaks.get(nmIx) + 1);
                            }
                        } else {
                            if (dist.firstPairDominates(weight, fbrWeight,
                                    weights.get(insertionIx), fbrWeights.get(insertionIx))) {
                                weights.set(insertionIx, weight);
                            }
                        }
                        for (int p = insertionIx - 1; p >= 0 &&
                                dist.firstPairStrictlyDominates(weight, fbrWeight, weights.get(p), fbrWeights.get(p)) &&
                                nmSigs.get(p).nmSum == currCount + currEpsilon;
                             --p) {
                            nmSigs.remove(p);
                            weights.removeAt(p);
                            fbrWeights.removeAt(p);
                            counts.remove(p);
                            nmStreaks.set(nmIx, nmStreaks.get(nmIx) - 1);
                            if(nmStreaks.get(nmIx) == 0){
                                nmStreaks.remove(nmIx);
                                nmRegions.remove(nmIx);
                            }
                        }
                    }
                }
                numElmts += prev.nmStreaks.get(i);
            }
        }
    }

    public static List<FbrSigCurve> sigCurves(Weights weights, FbrWeights fbrWeights,
                                              double minWeight, double minFbrWeight,
                                              int start, int end, Distances distance) {
        int k = end - start;
        Signature s = new Signature();
        List<FbrSigCurve> curves = new ArrayList<>(k);
        FbrSigCurve currCurve, prevCurve = new FbrSigCurve(minWeight, minFbrWeight);
        curves.add(prevCurve);
        for (int l = start; l < end; ++l) {
            currCurve = new FbrSigCurve(prevCurve, weights, fbrWeights, l, s, distance);
            curves.add(currCurve);
            prevCurve = currCurve;
        }

        return curves;
    }
    public int size() {
        return nmSigs.size();
    }
}


