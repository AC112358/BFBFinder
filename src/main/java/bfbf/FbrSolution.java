/*
 * This file is a part of the bfb java package for the analysis
 * of Breakage-Fusion-Bridge count vectors.
 *
 * Copyright (C) 2013 Shay Zakov, Marcus Kinsella, and Vineet Bafna.
 *
 * The bfb package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The bfb package is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact:
 * Shay Zakov:		zakovs@gmail.com
 */

package bfbf;

import bfbf.distances.Distances;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A data object class that stores a count vector, fbr vector, and some auxiliary meta data.
 *
 * @author Shay Zakov
 */
public class FbrSolution extends Solution {

    //protected int[] counts;
    protected int[] epsilons;
    //public Signature s;
    //protected double weight;
    protected double fbrWeight;
    //protected int totalCount;
    protected int nmSum;
    protected Distances dist;

    public FbrSolution() {
        counts = new int[1];
        epsilons = new int[1];
        s = new Signature();
        totalCount = 0;
        weight = 1;
        fbrWeight = 1;
        nmSum = -1;
        counts[0] = 0;
        epsilons[0] = 0;
        dist = new Distances();
    }


    public FbrSolution(int[] counts, int[] epsilons, Signature s, double weight, double fbrWeight) {
        super();
        this.counts = counts;
        this.epsilons = epsilons;
        totalCount = sum(counts);
        this.s = new Signature(s);
        this.weight = weight;
        this.fbrWeight = fbrWeight;
        nmSum = -1;
    }

    public FbrSolution(int[] counts, int[] epsilons, Signature s, double weight, double fbrWeight, Distances dist) {
        super();
        this.counts = counts;
        this.epsilons = epsilons;
        totalCount = sum(counts);
        this.s = new Signature(s);
        this.weight = weight;
        this.fbrWeight = fbrWeight;
        nmSum = -1;
        this.dist = dist;
    }

    private static int sum(int[] counts) {
        int totalCount = 0;
        for (int i = counts.length - 1; i >= 0; --i) {
            totalCount += counts[i];
        }
        return totalCount;
    }

    @Override
    public String toString() {
        String totalWeightString = "";
        if (dist != null){
            totalWeightString = " (combined weight: " + dist.combineWeights(weight, fbrWeight) + ")";
        }
        return ("{" + Arrays.toString(counts) + "," + Arrays.toString(getFbrCounts()) + "} (length: " + getLength() + "/" + counts.length + ", effective length: " +
                getEffectiveLength() + ", count weight: " + weight + ", fold-back weight: " + fbrWeight +
                totalWeightString + ")");
    }



    public int compareTo(FbrSolution other) {
        int nmSumDiff = compareToNm(other);
        if (nmSumDiff != 0){
            return nmSumDiff;
        }
        return s.lexCompare(other.s);
    }

    public int[] getFbrCounts(){
        int[] fbrCounts = new int[counts.length];
        for (int l = counts.length - 1; l >= 0; l--) {
            if (l == 0){
                fbrCounts[l] = counts[l] - 1 + epsilons[l];
                continue;
            }
            fbrCounts[l] = counts[l] - counts[l - 1] + epsilons[l] + epsilons[l - 1];
        }
        return fbrCounts;
    }

    public int[] displayFbrCounts(){
        int[] fbrCounts = getFbrCounts();
        for (int i = 0; i < fbrCounts.length; i++){
            fbrCounts[i] /= 2;
        }
        return fbrCounts;
    }

    public int compareToNm(FbrSolution other){
        if (nmSum == -1 && counts.length > 0){
            nmSum = counts[0] + epsilons[0];
        }
        if (other.nmSum == -1 && counts.length > 0){
            other.setNmSum(other.counts[0] + epsilons[0]);
        }
        int nmSumDiff = nmSum - other.nmSum;
        nmSum = -1;
        other.nmSum = -1;
        return nmSumDiff;
    }




    public double getFbrWeight() {
        return fbrWeight;
    }


    public void setFbrWeight(double fbrWeight) {
        this.fbrWeight = fbrWeight;
    }

    public int getNmSum(){
        if (nmSum == -1 && counts.length > 0){
            return counts[0] + epsilons[0];
        }
        return nmSum;
    }

    public void setNmSum(int nmSum){
        this.nmSum = nmSum;
    }




}
