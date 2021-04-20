package bfbf;

import bfbf.weights.Weights;
import bfbf.weights.FbrWeights;
import bfbf.distances.Distances;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.util.*;

public class FbrSigCurve{
    protected List<nmSignaturePair> nmSigs;
    //protected List<Signature> sigs;
    protected TDoubleList weights;
    protected TDoubleList fbrWeights;
    //protected List<Integer> nmSums;
    protected List<Integer> counts;
    protected List<Integer> epsilons;
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
        epsilons = new ArrayList<>(1);
        nmStreaks = new ArrayList<>(1);
        nmRegions = new ArrayList<>(1);
        weights.add(1);
        fbrWeights.add(1);
        counts.add(1);
        epsilons.add(1);
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

    /*private FbrSigCurve1(FbrSigCurve prev, Weights w, FbrWeights fw, int l, Signature s, Distances dist, int end) {
        nmSigs = new ArrayList<>(); // a series of lexicographically increasing signatures
        weights = new TDoubleArrayList(); // a series of decreasing weights
        fbrWeights = new TDoubleArrayList();
        counts = new ArrayList<>();
        epsilons = new ArrayList<>();
        this.minWeight = prev.minWeight;
        this.minFbrWeight = prev.minFbrWeight;
        nmStreaks = new ArrayList<>();
        nmRegions = new ArrayList<>();
        //boolean isValidFbr = (fw.getHeaviestCount(l) > -1);
        //boolean isPreValidFbr = (l < end - 1 && fw.getHeaviestCount(l + 1) > -1);
        boolean isValidFbr = (l > 0 && l < end - 1 && fw.getHeaviestCount(l) > -1);
        boolean isPreValidFbr = (l < end - 1 && fw.getHeaviestCount(l + 1) > -1);
        //System.out.println("l = " + l + " and isPreValidFbr == " + isPreValidFbr);
        //boolean isPostValidFbr = (l > 0 && fw.getHeaviestCount(l - 1) > -1);
        int maxCount = w.getMaxCount(l);
        for (int currCount = w.getMinCount(l); currCount <= maxCount; ++currCount) {
            //System.out.println("\t currCount = " + currCount);
            int numElmts = 0;
            double currWeight = w.getWeight(l, currCount);
            for (int i = 0; i < prev.nmStreaks.size(); i++) {
                int maxIx = prev.findSmallerWeightIx(minWeight / currWeight, numElmts, numElmts + prev.nmStreaks.get(i));
                int minIx = prev.findMaxSignature(currCount, numElmts, maxIx);
                for (int j = minIx; j < maxIx; ++j) {
                    s.setTo(prev.nmSigs.get(j).sig);
                    int prevNmSum = prev.nmSigs.get(j).nmSum;
                    int prevN = 0;
                    int prevM = 0;
                    prevN = prev.counts.get(j);
                    prevM = prev.epsilons.get(j);
                    //System.out.println("\t" + prev.nmSigs.get(j).sig + ": " + prevN + "," + prevM);
                    double weight = currWeight * prev.weights.get(j);
                    int minFbrCount = 0;
                    int maxFbrCount = currCount - 1;
                    if (isPreValidFbr) {
                        minFbrCount = 0;
                        maxFbrCount = currCount - 1;
                    }
                    if (isValidFbr) {
                        minFbrCount = fw.getMinCount(l, minFbrWeight, prevN, currCount, prevM);
                        maxFbrCount = fw.getMaxCount(l, minFbrWeight, prevN, currCount, prevM);
                        //System.out.println("\t minFbrCount = " + minFbrCount + ", maxFbrCount = " + maxFbrCount);
                    }
                    minFbrCount = Math.max(minFbrCount, 0);
                    if (l < end - 1){
                        maxFbrCount = Math.min(maxFbrCount, prevN - 1);
                    }
                    if (l == end - 1){
                        minFbrCount = currCount;
                        maxFbrCount = currCount;
                    }
                    //System.out.println("\t minFbrCount = " + minFbrCount + ", maxFbrCount = " + maxFbrCount);
                    for (int currEpsilon = minFbrCount; currEpsilon <= maxFbrCount; ++currEpsilon) {
                        //int tempFbrCount = currCount - prevN  + prevM + currEpsilon;
                        //tempFbrCount = tempFbrCount / 2;
                        if (tempFbrCount < minFbrCount || tempFbrCount > maxFbrCount){
                            //System.out.println("\t n = " + n +", m = " + m + ", tempFbrCount = " + tempFbrCount + ", maxFbrCount = " + maxFbrCount);
                            continue;
                        }
                        int tempEpsilon = currEpsilon;
                        s.setTo(prev.nmSigs.get(j).sig);
                        //System.out.println("\t currCount = " + currCount + ", currEpsilon = " + currEpsilon);
                        boolean toContinue = false;
                        if (s.minDecrement(currCount)){
                            //System.out.println("\t" + s);
                            if (!isValidFbr && !isPreValidFbr){
                                //tempEpsilon = s.numEpsilons;
                                //System.out.println(s + ": currCount = " + currCount + ", prevN = " + prevN +
                                //       ", numEpsilons = " + tempEpsilon);
                            }
                            if (!s.minDecrement(currCount - currEpsilon)){
                                toContinue = true;
                            }
                            //tempEpsilon += s.numEpsilons;
                        }else{
                            toContinue = true;
                        }
                        if (toContinue){
                            continue;
                        }

                        double currFbrWeight = 1.0;
                        if (isValidFbr) {
                            currFbrWeight = fw.getWeight(l, prevN, prevM, currCount, currEpsilon);
                        }
                        double fbrWeight = currFbrWeight * prev.fbrWeights.get(j);
                        int tempNmSum = currCount + tempEpsilon;
                        if (!isPreValidFbr){
                            tempNmSum = 0;
                        }
                        nmSignaturePair searchKey = new nmSignaturePair(tempNmSum, s);
                        //int insertionIx = Collections.binarySearch(sigs, s);
                        int nmIx = Collections.binarySearch(nmRegions, tempNmSum);
                        int insertionIx = Collections.binarySearch(nmSigs, searchKey);
                        //System.out.println(nmIx + "," + insertionIx);
                        if (nmIx < 0) {
                            insertionIx = -insertionIx - 1;
                            nmSigs.add(insertionIx, new nmSignaturePair(tempNmSum,
                                    new Signature(s)));
                            weights.insert(insertionIx, weight);
                            fbrWeights.insert(insertionIx, fbrWeight);
                            counts.add(insertionIx, currCount);
                            epsilons.add(insertionIx, tempEpsilon);
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
                                    !dist.firstPairDominates(weight, fbrWeight, weights.get(insertionIx), fbrWeights.get(insertionIx)
                                            )) {
                                // the next signature has a higher weight
                                continue;
                            } else {
                                nmSigs.add(insertionIx, new nmSignaturePair(tempNmSum,
                                        new Signature(s)));
                                weights.insert(insertionIx, weight);
                                fbrWeights.insert(insertionIx, fbrWeight);
                                counts.add(insertionIx, currCount);
                                epsilons.add(insertionIx, tempEpsilon);
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
                                nmSigs.get(p).nmSum ==  tempNmSum;
                             --p) {
                            nmSigs.remove(p);
                            weights.removeAt(p);
                            fbrWeights.removeAt(p);
                            counts.remove(p);
                            epsilons.remove(p);
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
*/
    private FbrSigCurve(FbrSigCurve prev, Weights w, FbrWeights fw, int l, Signature s, Distances dist, int end) {
        nmSigs = new ArrayList<>(); // a series of lexicographically increasing signatures
        weights = new TDoubleArrayList(); // a series of decreasing weights
        fbrWeights = new TDoubleArrayList();
        counts = new ArrayList<>();
        epsilons = new ArrayList<>();
        this.minWeight = prev.minWeight;
        this.minFbrWeight = prev.minFbrWeight;
        nmStreaks = new ArrayList<>();
        nmRegions = new ArrayList<>();
        //boolean isValidFbr = (fw.getHeaviestCount(l) > -1);
        //boolean isPreValidFbr = (l < end - 1 && fw.getHeaviestCount(l + 1) > -1);
        boolean isValidFbr = (l > 0 && l < end - 1 && fw.getHeaviestCount(l) > -1);
        boolean isPreValidFbr = (l < end - 1 && fw.getHeaviestCount(l + 1) > -1);
        //System.out.println("l = " + l + " and isPreValidFbr == " + isPreValidFbr);
        //boolean isPostValidFbr = (l > 0 && fw.getHeaviestCount(l - 1) > -1);
        int maxCount = w.getMaxCount(l);
        for (int currCount = w.getMinCount(l); currCount <= maxCount; ++currCount) {
            //System.out.println("currCount = " + currCount);
            int numElmts = 0;
            double currWeight = w.getWeight(l, currCount);
            for (int i = 0; i < prev.nmStreaks.size(); i++) {
                int maxIx = prev.findSmallerWeightIx(minWeight / currWeight, numElmts, numElmts + prev.nmStreaks.get(i));
                int minIx = prev.findMaxSignature(currCount, numElmts, maxIx);
                for (int j = minIx; j < maxIx; ++j) {
                    s.setTo(prev.nmSigs.get(j).sig);
                    int prevNmSum = prev.nmSigs.get(j).nmSum;
                    int prevN = 0;
                    int prevM = 0;
                    prevN = prev.counts.get(j);
                    prevM = prev.epsilons.get(j);
                    //System.out.println("\t" + prev.nmSigs.get(j).sig + ": " + prevN + "," + prevM);
                    double weight = currWeight * prev.weights.get(j);
                    int minFbrCount = 0;
                    int maxFbrCount = 0;
                    if (isPreValidFbr) {
                        minFbrCount = 0;
                        maxFbrCount = currCount - 1;
                    }
                    if (isValidFbr) {
                        minFbrCount = fw.getMinCount(l, minFbrWeight, prevN, currCount, prevM);
                        maxFbrCount = fw.getMaxCount(l, minFbrWeight, prevN, currCount, prevM);
                        //System.out.println("\t minFbrCount = " + minFbrCount + ", maxFbrCount = " + maxFbrCount);
                    }
                    minFbrCount = Math.max(minFbrCount, 0);
                    if (l < end - 1){
                        maxFbrCount = Math.min(maxFbrCount, currCount - 1);
                    }
                    if (l == end - 1){
                        minFbrCount = currCount;
                        maxFbrCount = currCount;
                    }
                    //System.out.println("\t minFbrCount = " + minFbrCount + ", maxFbrCount = " + maxFbrCount);
                    for (int currEpsilon = minFbrCount; currEpsilon <= maxFbrCount; ++currEpsilon) {
                        /*int tempFbrCount = currCount - prevN  + prevM + currEpsilon;
                        tempFbrCount = tempFbrCount / 2;
                        if (tempFbrCount < minFbrCount || tempFbrCount > maxFbrCount){
                            //System.out.println("\t n = " + n +", m = " + m + ", tempFbrCount = " + tempFbrCount + ", maxFbrCount = " + maxFbrCount);
                            continue;
                        }*/
                        int tempEpsilon = currEpsilon;
                        s.setTo(prev.nmSigs.get(j).sig);
                        //System.out.println("\t currCount = " + currCount + ", currEpsilon = " + currEpsilon);
                        boolean toContinue = false;
                        if (s.minDecrement(currCount)){
                            //System.out.println("\t" + s);
                            if (!isValidFbr && !isPreValidFbr){
                                tempEpsilon = s.numEpsilons;
                                //System.out.println(s + ": currCount = " + currCount + ", prevN = " + prevN +
                                //       ", numEpsilons = " + tempEpsilon);
                            }
                            if (!s.minDecrement(currCount - currEpsilon)){
                                toContinue = true;
                            }
                            tempEpsilon += s.numEpsilons;
                        }else{
                            toContinue = true;
                        }
                        if (toContinue){
                            continue;
                        }

                        double currFbrWeight = 1.0;
                        if (isValidFbr) {
                            currFbrWeight = fw.getWeight(l, prevN, currCount, prevM, currEpsilon);
                            //System.out.println("\t currFbrWeight = " + currFbrWeight);
                            //System.out.println("\t currEpsilon = " + currEpsilon);
                        }
                        double fbrWeight = currFbrWeight * prev.fbrWeights.get(j);
                        int tempNmSum = currCount + tempEpsilon;
                        if (!isPreValidFbr){
                            tempNmSum = 0;
                        }

                        nmSignaturePair searchKey = new nmSignaturePair(tempNmSum, s);
                        //int insertionIx = Collections.binarySearch(sigs, s);
                        int nmIx = Collections.binarySearch(nmRegions, tempNmSum);
                        int insertionIx = Collections.binarySearch(nmSigs, searchKey);
                        //System.out.println(nmIx + "," + insertionIx);
                        if (nmIx < 0) {
                            insertionIx = -insertionIx - 1;
                            nmSigs.add(insertionIx, new nmSignaturePair(tempNmSum,
                                    new Signature(s)));
                            weights.insert(insertionIx, weight);
                            fbrWeights.insert(insertionIx, fbrWeight);
                            counts.add(insertionIx, currCount);
                            epsilons.add(insertionIx, tempEpsilon);
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
                                    !dist.firstPairDominates(weight, fbrWeight, weights.get(insertionIx), fbrWeights.get(insertionIx)
                                    )) {
                                // the next signature has a higher weight
                                continue;
                            } else {
                                nmSigs.add(insertionIx, new nmSignaturePair(tempNmSum,
                                        new Signature(s)));
                                weights.insert(insertionIx, weight);
                                fbrWeights.insert(insertionIx, fbrWeight);
                                counts.add(insertionIx, currCount);
                                epsilons.add(insertionIx, tempEpsilon);
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
                                nmSigs.get(p).nmSum ==  tempNmSum;
                             --p) {
                            nmSigs.remove(p);
                            weights.removeAt(p);
                            fbrWeights.removeAt(p);
                            counts.remove(p);
                            epsilons.remove(p);
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
            currCurve = new FbrSigCurve(prevCurve, weights, fbrWeights, l, s, distance, end);
            curves.add(currCurve);
            prevCurve = currCurve;
            //System.out.println(currCurve.nmSigs);
            //System.out.println(currCurve.counts);
            //System.out.println(currCurve.epsilons);
        }
        return curves;
    }
    public int size() {
        return nmSigs.size();
    }
}


