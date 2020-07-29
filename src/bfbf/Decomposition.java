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

import bfbf.palindromes.*;

/**
 * An object representing a decomposition of an l-BFB palindrome collection.
 *
 * @author Shay Zakov
 */
public class Decomposition {

    protected int size, r;
    PalindromeCollection[] L, H;
    int l;

    public Decomposition(int maxR, int l) {
        size = 0;
        r = 0;
        L = new PalindromeCollection[maxR + 1];
        H = new PalindromeCollection[maxR + 1];
        this.l = l;

        for (int i = 0; i <= maxR; ++i) {
            L[i] = new PalindromeCollection();
            H[i] = new PalindromeCollection();
        }
    }

    public int[] signature() {
        int[] s = new int[r + 2];

        s[0] = L[0].fullSize();
        for (int i = 1; i <= r; ++i) {
            s[i] = L[i].fullSize() - L[i - 1].fullSize() - H[i - 1].fullSize() / 2 + Math.max(s[i - 1], 0);
        }

        return s;
    }

    public void wrap() {
        for (int i = 0; i < r; ++i) {
            L[i].wrap();
            H[i].wrap();
        }
        --l;
    }

    public void addEmpty(int m) {
        if (m > 0) {
            int d = Env.dig(m);
//			EmptyPalindrome epsilon = new EmptyPalindrome(l);
            L[d].add(EmptyPalindrome.SINGLETON);
            int factor = 2;
            for (int i = d + 1; i < r; ++i, factor *= 2) {
                for (int j = 0; j < H[i].uniqueSize(); ++j) {
                    H[d].add(H[i].get(j), factor * H[i].getMultiplicity(j));
                }
                H[i].clear();

                for (int j = 0; j < L[i].uniqueSize(); ++j) {
                    H[d].add(L[i].get(j), factor * L[i].getMultiplicity(j));
                }
                L[i].clear();
            }
            H[d].add(EmptyPalindrome.SINGLETON, (int) (m / Math.pow(2, d) - 1));

            r = d + 1;
            size += m;
        }
    }

    public void rightFold() {
        int g = r - 1;
        size -= Math.pow(2, r);

        ConvexPalindrome gamma = ConvexPalindrome.make((BFBPalindrome) L[g].removeMin(), null);

        // Accumulating the convexed l-palindrome:
        for (--g; H[g].isEmpty(); --g) {
            gamma = ConvexPalindrome.make((BFBPalindrome) L[g].removeMin(), gamma);
        }

        for (; H[r - 1].isEmpty() && L[r - 1].isEmpty(); --r) ; // updating r(B)

        gamma = ConvexPalindrome.make((BFBPalindrome) L[g].removeMin(), gamma);

        Block beta = (Block) H[g].removeMin();
        H[g].removeMin(); // need to remove two copies of beta from H[g].
//		Collections.reverse(gamma.palList);
        L[g].add(CompositPalindrome.make(gamma, beta));
    }

    public void fold(int n) {
        if (n != size) {

            if (n > size) {
                addEmpty(n - size);
            } else {
                int dMax = Env.dig(n - size), dOpt = 0;
                int s_d = L[0].fullSize();
                int Delta_d = 0;

                for (int d = 0, factor = 1; d <= dMax; ++d, factor *= 2) {
                    // Invariant: factor = 2^d, s_d is the d-th element in the
                    // signature of the collection, and Delta_d is the d-th element
                    // in the corresponding series Delta.

                    if (n >= factor * (Math.max(s_d + 1, 0)) + Delta_d) {
                        dOpt = d;
                    }
                    // Setting Delta_d and s_d with respect to d+1:
                    Delta_d += factor * (Math.abs(s_d));
                    s_d = L[d + 1].fullSize() - L[d].fullSize() - H[d].fullSize() / 2 + Math.max(s_d, 0);
                }


                addEmpty((int) Math.pow(2, dOpt));

                while (size > n) {
                    rightFold();
                }

                addEmpty(n - size);
            }
        }
    }

}
