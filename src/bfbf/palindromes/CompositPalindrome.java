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

import java.util.ArrayList;
import java.util.List;


/**
 * An composite implementation of the abstract Palindrome class.
 *
 * @author Shay Zakov
 */
public class CompositPalindrome extends BFBPalindrome {

    private static final CompositPalindrome cashKey = new CompositPalindrome(null, null);
    private static final List<Palindrome> pool = new ArrayList<Palindrome>();
    private static final PalindromeFactory<CompositPalindrome> factory =
            new CompositeFactory();


    protected CompositPalindrome(ConvexPalindrome convex, Block beta) {
        super(convex, beta);
    }

    public static CompositPalindrome make(ConvexPalindrome convex, Block beta) {
        return make(convex, beta, cashKey, factory);
    }

    @Override
    protected int classCompareTo(Palindrome other) {
        if (other instanceof CompositPalindrome) return super.classCompareTo(other);
        else return 1;
    }

    public CompositPalindrome extendCenter(BFBPalindrome toAdd) {
        if (toAdd.depth() > depth()) {
            throw new IllegalArgumentException("Center depth must be " +
                    "lower than or equals to wrap depth.");
        }

        ConvexPalindrome convex = ConvexPalindrome.make(toAdd, null);
        if (isEmpty(center)) {
            return make(convex, (Block) wrap);
        } else {
            return make(((ConvexPalindrome) center).extend(convex), (Block) wrap);
        }
    }

    public int nestingDeg() {
        if (center == null) {
            return 0;
        } else return ((ConvexPalindrome) center).nestingDeg();
    }

    @Override
    public int centerDeg() {
        if (isEmpty(center) || center.depth() < depth()) {
            return 1;
        } else {
            return 1 + center.centerDeg();
        }
    }

    @Override
    public int minPartDepth() {
        if (isEmpty(center)) return wrap.depth();
        else return center.minPartDepth();
    }


    @Override
    public List<Palindrome> pool() {
        return pool;
    }


    @Override
    protected Palindrome primary() {
        return wrap;
    }

    @Override
    protected Palindrome secondary() {
        return center;
    }


    public static class CompositeFactory extends PalindromeFactory<CompositPalindrome> {

        @Override
        public CompositPalindrome make(Palindrome center, Palindrome wrap) {
            return new CompositPalindrome((ConvexPalindrome) center, (Block) wrap);
        }

    }
}
