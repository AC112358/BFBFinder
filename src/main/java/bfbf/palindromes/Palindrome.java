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

package bfbf.palindromes;

import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * An abstract palindrome class.
 * <p>
 * A palindrome is defined by a center palindrome P_c and a wrap
 * palindrome P_w, and represents the structure P_w P_c P_w.
 *
 * @author Shay Zakov
 */
public abstract class Palindrome implements Comparable<Palindrome> {

    //	protected int depth;  // The number of nested layers.

    protected static final Map<Palindrome, Palindrome> cashedInstances =
            new HashMap<>();
    protected static int cashedDepth = 5; // palindromes of depth at most cashedDepth will be cashed for re-usage.
    protected Palindrome center, wrap;

    /**
     * Constractor of a palindrome based on its center and wrap components.
     *
     * @param center the center palindrome
     * @param wrap   the wrap palindrome
     */
    protected Palindrome(Palindrome center, Palindrome wrap) {
        this.center = center;
        this.wrap = wrap;
    }

    protected static <P extends Palindrome> P make(Palindrome center, Palindrome wrap,
                                                   P cashKey, PalindromeFactory<P> factory) {
        cashKey.center = center;
        cashKey.wrap = wrap;
        Palindrome p = cashedInstances.get(cashKey);
        if (p == null) {
            List<Palindrome> pool = cashKey.pool();
            if (!pool.isEmpty()) {
                p = pool.remove(pool.size() - 1);
                p.center = center;
                p.wrap = center;
            } else {
                p = factory.make(center, wrap);
            }
            if (p.depth() <= cashedDepth) {
                cashedInstances.put(p, p);
            }
        }
        return (P) p;
    }

    /**
     * Computes the difference between depths of two (possibly null) palindromes.
     * Note that the depth of "null" is defined to be -1, whereas the depth of an
     * empty palindrome is 0.
     *
     * @param first  First palindrome
     * @param second Second palindrome.
     * @return d1-d2, where d1 os the depth of the first palindrome and d2 is the
     * depth of the second palindrome.
     */
    protected static int depthDiff(Palindrome first, Palindrome second) {
        int firstDepth, secondDepth;

        if (first == null) firstDepth = -1;
        else firstDepth = first.depth();
        if (second == null) secondDepth = -1;
        else secondDepth = second.depth();

        return firstDepth - secondDepth;
    }

    public static boolean isEmpty(Palindrome p) {
        return p == null || p.length() == 0;
    }

    protected static void retain(Palindrome p) {
        if (!isEmpty(p) && !cashedInstances.containsKey(p)) {
            retain(p.center);
            retain(p.wrap);
            p.center = null;
            p.wrap = null;
            p.pool().add(p);
        }
    }

    /**
     * Returns the palindrome's depth, i.e. the number of unique characters
     * in the string represented by this palindrome.
     *
     * @return the depth of the palindrome
     */
    public int depth() {
        Palindrome primary = primary();
        if (primary == null) return 0;
        else return primary.depth();
    }

    /**
     * Constructs the sequence segment indices represented by this
     * palindrome.
     *
     * @return the palindrome's segment sequence
     */
    public int[] seq() {
        return seq(0);
    }

    /**
     * Constructs the sequence of segment indices represented by this
     * palindrome, given a minimal (base) segment index.
     *
     * @param base the index of the minimal segment
     * @return the palindrome's segment sequence
     */
    public int[] seq(int base) {
        int[] seq = new int[length()];
        fillSeq(seq, base, 0);
        return seq;
    }

    protected int fillSeq(int[] seq, int base, int startIx) {
        int length = length();
        int wrapLength = 0;
        if (wrap != null) {
            wrapLength = wrap.fillSeq(seq, base, startIx);
            System.arraycopy(seq, startIx, seq, startIx + length - wrapLength, wrapLength);
        }
        if (center != null) {
            center.fillSeq(seq, base, startIx + wrapLength);
        }
        return length;
    }

    @Override
    public int compareTo(Palindrome other) {
        if (other == this) {
            return 0;
        } else if (other == null) {
            return -1;
        } else if (other.depth() != depth()) {
            return other.depth() - depth();
        } else {
            return classCompareTo(other);
        }
    }

    // TODO: possibly inconsistent comparisons between convex and composite
    //  palindromes (return value depends on the order).
    protected int classCompareTo(Palindrome other) {
        int res = primary().compareTo(other.primary());
        if (res == 0) {
            if (secondary() != null) {
                if (other.secondary() == null) return 1;
                else if (secondary() != primary() || other.secondary() != other.primary()) {
                    res = secondary().compareTo(other.secondary());
                }
            } else if (other.secondary() != null) {
                res = -1;
            }
        }
        return res;
    }

    public String halfString() {
        String fullString = toString();
        return fullString.substring(0, fullString.length() / 2);
    }

    @Override
    public String toString() {
        return toString('A');
    }

    public String toString(char base) {
        int[] seq = seq();
        if (seq.length == 0) {
            return "e";
        }

        char[] chars = new char[seq.length];
        for (int i = 0; i < seq.length; ++i) {
            chars[i] = (char) (base + seq[i]);
        }
        return new String(chars);
    }

    /**
     * This method computes the "center degree" of a palindrome, which is
     * used for enumerating all foldings of a BFB palindrome collection.
     * <p>
     * The center degree of a convex palindrome is defined to be the center
     * degree of its (BFB palindrome) center.
     * The center degree of a block is defined to be zero. For a composite
     * palindrome, the center degree is zero if its center is strictly lower
     * than its wrap, and otherwise it is 1 plus the center degree of its
     * center.
     *
     * @return The center degree of the BFB palindrome.
     */
    public abstract int centerDeg();

    // TODO: document!
    public abstract int minPartDepth();

    /**
     * The length of the segment sequence corresponding to this palindrome.
     *
     * @return the length of the segment sequence corresponding to this palindrome.
     */
    public int length() {
        int length = 0;
        if (center != null) {
            length = center.length();
        }
        if (wrap != null) {
            length += 2 * wrap.length();
        }
        return length;
    }

    public abstract List<Palindrome> pool();


    /**
     * @return the deeper palindrome among this palindrome's center and wrap
     */
    protected abstract Palindrome primary();

    /**
     * @return the shallower palindrome among this palindrome's center and wrap
     */
    protected abstract Palindrome secondary();

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((center == null) ? 0 : center.hashCode());
        result = prime * result + ((wrap == null) ? 0 : wrap.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Palindrome other = (Palindrome) obj;
        if (center == null) {
            if (other.center != null)
                return false;
        } else if (!center.equals(other.center))
            return false;
        if (wrap == null) {
            if (other.wrap != null)
                return false;
        } else if (!wrap.equals(other.wrap))
            return false;
        return true;
    }

    /**
     * @return the segment count vector corresponding to this palindrome
     */
    public int[] counts() {
        int[] counts = new int[depth()];
        fillCounts(counts, 0, 1);
        return counts;
    }

    protected void fillCounts(int[] counts, int base, int factor) {
        if (center != null) {
            center.fillCounts(counts, base, factor);
        }
        if (wrap != null) {
            wrap.fillCounts(counts, base, 2 * factor);
        }
    }

    public abstract static class PalindromeFactory<P extends Palindrome> {

        public abstract P make(Palindrome center, Palindrome wrap);

    }
}
