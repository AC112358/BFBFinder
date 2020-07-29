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

package bfbf.instance_generation;

import java.util.Arrays;
import java.util.Iterator;

/**
 * An iterator over all count vectors of certain length k and maximum entry
 * value n.
 *
 * @author Shay Zakov
 */
public class CVIterator implements Iterator<int[]> {

    private int[] next;
    private int k, n, currSum, maxFirst;

    public CVIterator(int k, int n) {
        assert n >= k;
        this.k = k;
        this.n = n;
        maxFirst = n - k + 1;
        next = new int[k];
        Arrays.fill(next, 1);
        currSum = k - 2;
        next[k - 2] = 0;
    }

    @Override
    public boolean hasNext() {
        return next[0] < maxFirst;
    }

    @Override
    public int[] next() {
        int l = k - 2;

        if (currSum == n - 1) {
            for (; next[l] == 1; --l) {
                next[l] = 1;
            }
            currSum -= next[l] - 1;
            next[l] = 1;
            --l;
        }

        ++next[l];
        ++currSum;
        next[k - 1] = n - currSum;
        return next;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
