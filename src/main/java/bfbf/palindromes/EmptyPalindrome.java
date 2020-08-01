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

import java.util.List;

/**
 * An empty palindrome implementation of the abstract Palindrome class.
 *
 * @author Shay Zakov
 */
public class EmptyPalindrome extends BFBPalindrome {

    public static final EmptyPalindrome SINGLETON = new EmptyPalindrome();

    private EmptyPalindrome() {
        super(null, null);
    }

    public boolean equals(Object other) {
        return other == SINGLETON;
    }

    @Override
    public int centerDeg() {
        return 0;
    }

    @Override
    public int minPartDepth() {
        return 0;
    }

    @Override
    public List<Palindrome> pool() {
        return null;
    }

    @Override
    protected Palindrome primary() {
        return null;
    }

    @Override
    protected Palindrome secondary() {
        return null;
    }


}
