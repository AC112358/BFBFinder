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
public class Solution implements Comparable<Solution> {

    public int[] s;
    public double error;
    public int totalCount;
    int[] counts;

    public Solution() {
    }


    public Solution(int[] counts, int[] s, double error) {
        super();
        this.counts = counts;
        totalCount = sum(counts);
        this.s = s;
        this.error = error;
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
        if (counts == null || counts.length == 0) {
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

    public static int start(int[] vec) {
        int start = 0;
        for (; start < vec.length && vec[start] == 0; ++start) ;
        return start;
    }

    public static int end(int[] vec) {
        int end = vec.length;
        for (; end > 0 && vec[end - 1] == 0; --end) ;
        return end;
    }

    @Override
    public int compareTo(Solution other) {
        int i = 0;
        int last = Math.min(s.length, other.s.length);
        for (; i < last && s[i] == other.s[i]; ++i) ;
        if (i == last) {
            return Env.errorModel.compareErrors(error, other.error);
        } else return s[i] - other.s[i];
    }

    public int getS(int i) {
        if (i < s.length) return s[i];
        else return 0;
    }

    @Override
    public String toString() {
        return Arrays.toString(counts) + " (length: " + getLength() + "/" + counts.length + ", effective length: " + getEffectiveLength() + ", error: " + error + ")";// Env.errorModel.normlizeError(error) +")";
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
        int r = s.length - 1;
        for (; r > 0 && s[r - 1] == 0; --r) ;
        if (r > 0 && s[r - 1] < 0) {
            --r;
        }
        return r;
    }

    public int getCount(int i) {
        return counts[i];
    }

    public int[] getEffectiveCounts() {
        return effectiveCounts(counts);
    }

}
