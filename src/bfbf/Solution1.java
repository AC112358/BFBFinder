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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A data object class that stores a count vector and some auxiliary meta data.
 *
 * @author Shay Zakov
 */
public class Solution1 implements Comparable<Solution1> {

    protected int[] counts;
    protected Signature s;
    protected double weight;
    protected int totalCount;

    public Solution1() {
        counts = new int[0];
        s = new Signature();
        totalCount = 0;
        weight = 1;
    }


    public Solution1(int[] counts, Signature s, double weight) {
        super();
        this.counts = counts;
        totalCount = sum(counts);
        this.s = new Signature(s);
        this.weight = weight;
    }

    private static int sum(int[] counts) {
        int totalCount = 0;
        for (int i = counts.length - 1; i >= 0; --i) {
            totalCount += counts[i];
        }
        return totalCount;
    }

    public static int[] effectiveCounts(int[] counts) {
        int prev = 0;
        List<Integer> ls = new ArrayList<>();
        for (int i = 0; i < counts.length; ++i) {
            if (counts[i] != prev) {
                ls.add(counts[i]);
                prev = counts[i];
            }
        }
        if (counts[counts.length - 1] == 0) {
            ls.remove(ls.size() - 1);
        }

        int[] effectiveCounts = new int[ls.size()];
        for (int j = 0; j < ls.size(); ++j) {
            effectiveCounts[j] = ls.get(j);
        }

        return effectiveCounts;
    }

    public static int effectiveLength(int[] counts) {
        if (counts.length == 0) {
            return 0;
        }
        int length = 0, prev = 0;
        for (int i = counts.length - 1; i >= 0; --i) {
            if (counts[i] != prev) {
                ++length;
                prev = counts[i];
            }
        }
        if (counts[0] != 0) {
            return length;
        } else {
            return length - 1;
        }
    }

    @Override
    public int compareTo(Solution1 other) {
        int res = s.compareTo(other.s);
//		TODO: do we also want to compare weights???
        if (res == 0) {
            return (int) Math.signum(other.weight - weight);
        }
        return res;
    }

    public int getS(int i) {
        return s.get(i);
    }

    @Override
    public String toString() {
        return Arrays.toString(counts) + " (length: " + getLength() + "/" + counts.length + ", effective length: " + getEffectiveLength() + ", weight: " + weight + ")";// Env.errorModel.normlizeError(error) +")";
    }

    /**
     * @return the length of the vector after trimming zeros at both ends.
     */
    public int getLength() {
        int i = 0, j = counts.length;
        for (; i < j && counts[i] == 0; ++i) ;
        for (; j > i && counts[j - 1] == 0; --j) ;
        return j - i;
    }

    /**
     * @return the length of the vector after trimming zeros at both ends and
     * regarding consecutive identical counts as one count.
     */
    public int getEffectiveLength() {
        return effectiveLength(counts);
    }

    public int getR() {
        return s.r();
    }

    public int getCount(int i) {
        return counts[i];
    }

    public int[] getEffectiveCounts() {
        return effectiveCounts(counts);
    }


    public double getWeight() {
        return weight;
    }


    public void setWeight(double weight) {
        this.weight = weight;
    }


}
