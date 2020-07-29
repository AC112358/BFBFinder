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
 * An l-block implementation of the abstract Palindrome class.
 *
 * @author Shay Zakov
 */
public class Block extends BFBPalindrome {

    private static final Block cashKey = new Block(null);
    private static final List<Palindrome> pool = new ArrayList<Palindrome>();
    private static final PalindromeFactory<Block> factory = new BlockFactory();

    protected Block(Palindrome internal) {
        super(internal, null);
    }

    public static Block make(BFBPalindrome center) {
        return make(center, null, cashKey, factory);
    }

    @Override
    protected void fillCounts(int[] counts, int base, int factor) {
        counts[base] += factor;
        if (center != null) {
            center.fillCounts(counts, base + 1, factor);
        }
    }

    @Override
    protected int classCompareTo(Palindrome other) {
        if (other instanceof Block) return center.compareTo(other.center);
        else return -1;
    }

    @Override
    protected int fillSeq(int[] seq, int base, int startIx) {
        int length = length();
        seq[startIx] = base;
        seq[startIx + length - 1] = base;
        center.fillSeq(seq, base + 1, startIx + 1);
        return length;
    }

    @Override
    public int centerDeg() {
        return 0;
    }

    @Override
    public int minPartDepth() {
        return depth();
    }

    @Override
    public int length() {
        return center.length() + 2;
    }

    @Override
    public int depth() {
        return center.depth() + 1;
    }

    @Override
    public List<Palindrome> pool() {
        return pool;
    }

    @Override
    protected Palindrome primary() {
        return center;
    }

    @Override
    protected Palindrome secondary() {
        return center;
    }

    public static class BlockFactory extends PalindromeFactory<Block> {

        @Override
        public Block make(Palindrome center, Palindrome wrap) {
            return new Block(center);
        }

    }
}
