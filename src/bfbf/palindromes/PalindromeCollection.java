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

import bfbf.*;
import bfbf.weights.NoErrorModel;
import bfbf.weights.PoissonErrorModel;
import bfbf.weights.Weights;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


/**
 * A collection of l-palindromes, sorted with decreasing depths.
 *
 * @author Shay Zakov
 */
public class PalindromeCollection extends ArrayList<BFBPalindrome> {

    /**
     *
     */
    private static final long serialVersionUID = 7627792656453386450L;

    private final List<Integer> multiplicities;
    private int fullSize;

    public PalindromeCollection() {
        multiplicities = new ArrayList<>();
        fullSize = 0;
    }


    public PalindromeCollection(PalindromeCollection toCopy) {
        super(toCopy);
        multiplicities = new ArrayList<>(toCopy.multiplicities);
        fullSize = toCopy.fullSize;
    }

    public static void main(String[] args) {

        test508();
        System.exit(0);

        String[] strings = {
                "CBBBBCDEEDCBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABCDE",
                "CBAABBAABCDEEDCBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABCDE",
                "CBAAAABBAAAABCDEEDCBAAAAAAAAAAAAAAAAAAAAAAAAAABCDE",
                "CBAAAAAABBAAAAAABCDEEDCBAAAAAAAAAAAAAAAAAAAAAABCDE",
                "CBAAAAAAAABBAAAAAAAABCDEEDCBAAAAAAAAAAAAAAAAAABCDE",
                "CBAABCDEEDCBAAAAAAAAAAAAAAAABBAAAAAAAAAAAAAAAABCDE",
                "CCDEEDCBAAAAAAAAAAAAAAAABBAABBAAAAAAAAAAAAAAAABCDE",
                "CBAAAAAABCDEEDCBAAAAAAAAAAAAAABBAAAAAAAAAAAAAABCDE",
                "CCDEEDCBAAAAAAAAAAAAAABBAAAAAABBAAAAAAAAAAAAAABCDE",
                "CBAAAAAAAAAABBAAAAAAAAAABCDEEDCBAAAAAAAAAAAAAABCDE",
                "CBAAAAAAAAAABCDEEDCBAAAAAAAAAAAABBAAAAAAAAAAAABCDE",
                "CCDEEDCBAAAAAAAAAAAABBAAAAAAAAAABBAAAAAAAAAAAABCDE",
                "CCDEEDCBAAAAAAAAAABBAAAAAAAAAAAAAABBAAAAAAAAAABCDE",
                "CBAAAAAAAAAAAAAABCDEEDCBAAAAAAAAAABBAAAAAAAAAABCDE",
                "CBAAAAAAAAAAAABBAAAAAAAAAAAABCDEEDCBAAAAAAAAAABCDE",
                "CCDEEDCBAAAAAAAABBAAAAAAAAAAAAAAAAAABBAAAAAAAABCDE",
                "CBAAAAAAAAAAAAAAAAAABCDEEDCBAAAAAAAABBAAAAAAAABCDE",
                "CCDEEDCBAAAAAABBAAAAAAAAAAAAAAAAAAAAAABBAAAAAABCDE",
                "CBAAAAAAAAAAAAAAAAAAAAAABCDEEDCBAAAAAABBAAAAAABCDE",
                "CBAAAAAAAAAAAAAABBAAAAAAAAAAAAAABCDEEDCBAAAAAABCDE",
                "CCDEEDCBAAAABBAAAAAAAAAAAAAAAAAAAAAAAAAABBAAAABCDE",
                "CBAAAAAAAAAAAAAAAAAAAAAAAAAABCDEEDCBAAAABBAAAABCDE",
                "CCDEEDCBAABBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBAABCDE",
                "CBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABCDEEDCBAABBAABCDE",
                "CBAAAAAAAAAAAAAAAABBAAAAAAAAAAAAAAAABCDEEDCBAABCDE"
        };

        String[] original = reverse(strings, 'A', 'E');

        // ABBAAAABBA
        // ABBAABBAAA
        // ABBBBAA
        // ABBAABB
        //				int[] counts = {3, 3, 4, 6, 34}; //34,6,4,3,3
        int[] counts = {5, 3, 14, 6, 12, 6, 10, 8, 12, 8};
        //		int[] counts = {13, 11, 14, 10, 12, 14, 10, 14, 4}; /// this vector has 5,797,288 BFB strings //{13, 11, 14, 10, 12, 14, 10, 14, 4};


        Weights w = new NoErrorModel().getWeights(counts, 1);
        //		List<BFBPalindrome> palindromes =
        BFB.allBFBStrings(w, 1, counts.length);
        //		System.out.println(Bfb.countAllBFBStrings(w, 1, counts.length, 2005));
        BFB.allBFBPairwiseProjections(w, 1, counts.length);

        //		System.out.println(palindromes.size());
        System.exit(0);
        //		String prevString = palindromes.get(0).halfString();
        //
        //		for (int i=1; i<palindromes.size(); ++i){
        //			String currString = palindromes.get(i).halfString();
        //			System.out.println(currString);
        //			System.out.println(Arrays.equals(counts, Env.count(currString)) + ", "
        //					+ Bfb.isBFB(currString));
        //			if (!Arrays.equals(counts, Env.count(currString)) ||
        //					currString.equals(prevString) || !Bfb.isBFB(currString)){
        //				System.err.println("BUG");
        //			}
        //			prevString = currString;
        //			//			System.out.println(original[i]);
        //			//			System.out.println(solStrings[i].equals(original[i]));
        //			System.out.println();
        //		}
    }

    private static String[] reverse(String[] strings, char bottom, char top) {
        String[] res = new String[strings.length];
        for (int i = 0; i < strings.length; ++i) {
            String reversed = "";
            for (int j = strings[i].length() - 1; j >= 0; --j) {
                reversed += (char) (bottom + top - strings[i].charAt(j));
            }
            res[i] = reversed;
        }
        return res;
    }

    public static void test508() {
        double minWeight = 0.8;
        int[] counts = {5, 3, 14, 5, 12, 5, 11, 7, 12, 7, 14, 4, 14, 4, 14, 5, 3, 5, 3};
        PoissonErrorModel em = new PoissonErrorModel();
        Weights w = em.getWeights(counts, minWeight);
        int minLength = 9;

        System.out.println("Input: " + Arrays.toString(counts));
        System.out.println("Number of input segments: " + counts.length);
        System.out.println("Error model: " + em.toString() + ", minimum solution weight: " + minWeight);

        List<int[]> allHeavyBFBSubVectors = Signature.allHeavyBFBSubVectors(w, minWeight, minLength);
        System.out.println("Number of heavy BFB sub-vectors: " + allHeavyBFBSubVectors.size());
        System.out.println();

        int stringSolutions = 0;
        NoErrorModel nem = new NoErrorModel();

        for (int[] vec : allHeavyBFBSubVectors) {
            System.out.println(Arrays.toString(vec));
            vec = Solution1.effectiveCounts(vec);
            w = nem.getWeights(vec, 1);
            //			List<BFBPalindrome> palindromes =
            //			Bfb.allBFBStrings(w, 1, counts.length);
            int x = BFB.countAllBFBStrings(w, 1, vec.length);
            System.out.println("Number of BFB strings: " + x + "\n");
            stringSolutions += x;
            //			Bfb.allBFBPairwiseProjections(w,1,counts.length);
        }

        System.out.println("Total number of BFB strings: " + stringSolutions + "\n");
    }

    @Override
    public boolean add(BFBPalindrome p) {
        return add(p, 1);
    }

    public boolean add(BFBPalindrome p, int multiplicity) {
        if (multiplicity > 0) {
            int ix = Collections.binarySearch(this, p);
            if (ix < 0) {
                ix = -ix - 1;
                super.add(ix, p);
                multiplicities.add(ix, multiplicity);
            } else {
                multiplicities.set(ix, multiplicities.get(ix) + multiplicity);
            }
            fullSize += multiplicity;
        }
        return true;
    }

    public int uniqueSize() {
        return super.size();
    }

    /**
     * Wrapping all elements in the collection (rendering it into a
     * collection of blocks).
     */
    public void wrap() {
        for (int i = uniqueSize() - 1; i >= 0; --i) {
            set(i, get(i).wrap());
        }
    }

    /**
     * The inverse of the "wrap" method - given this is a collection
     * of blocks, replacing each block in the collection by its center
     * BFB palindrome.
     */
    public void unwrap() {
        for (int i = uniqueSize() - 1; i >= 0; --i) {
            set(i, (BFBPalindrome) get(i).center);
        }
    }

    @Override
    public PalindromeCollection clone() {
//		PalindromeCollection clone = new PalindromeCollection();
//		clone.addAll(this);
        PalindromeCollection clone = (PalindromeCollection) super.clone();
        clone.multiplicities.addAll(this.multiplicities);
        clone.fullSize = this.fullSize;
        return clone;
    }

    // TODO: no usages - consider removing
    public int totalLength() {
        int length = 0;
        for (int i = uniqueSize() - 1; i >= 0; --i) {
            length += multiplicities.get(i) * get(i).length();
        }
        return length;
    }

    /**
     * @return the total segment counts in all members of the collection
     */
    public int[] counts() {
        int k = depth();
        int[] counts = new int[k];
        for (int i = uniqueSize() - 1; i >= 0; --i) {
            Palindrome p = get(i);
            p.fillCounts(counts, 0, getMultiplicity(i));
        }
        return counts;
    }

    /**
     * @return the max depth of an element in the collection,
     * or 0 if the collection is empty.
     */
    public int depth() {
        if (!isEmpty()) {
            return get(0).depth();
        } else return 0;
    }

    public int getMultiplicity(int i) {
        return multiplicities.get(i);
    }

    /**
     * Multiples the number of occurences of each element in the collection by
     * the given factor.
     *
     * @param factor the multiplication factor
     */
    public void multiply(double factor) {
        for (int i = multiplicities.size() - 1; i >= 0; --i) {
            multiplicities.set(i, (int) (multiplicities.get(i) * factor));
        }
    }

    /**
     * Removes all elements from the collection.
     */
    public void clear() {
        super.clear();
        multiplicities.clear();
        fullSize = 0;
    }

    /**
     * Removes an element with minimum depth from the collection.
     *
     * @return the removed element
     */
    public Palindrome removeMin() {
        Palindrome min;
        int last = uniqueSize() - 1;
        int minMultiplicity = multiplicities.get(last);
        if (minMultiplicity == 1) {
            min = remove(last);
            multiplicities.remove(last);
        } else {
            min = get(last);
            multiplicities.set(last, minMultiplicity - 1);
        }
        --fullSize;
        return min;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("{");
        if (uniqueSize() > 0) {
            sb.append(multiplicities.get(0)).append(" * ").append(get(0));
            for (int i = 1; i < uniqueSize(); ++i) {
                sb.append(", ").append(multiplicities.get(i)).append(" * ").append(get(i));
            }
        }
        sb.append("}");

        return sb.toString();
    }


    //	-------------------------------------------------------------------

    /**
     * Exhaustively generating all possible foldings of the current collection,
     * which weights and signatures are valid with respect to a given signature
     * oracle.
     *
     * @param foldings       TODO
     * @param foldingWeights TODO
     * @param w              a weight function.
     * @param l              the index of the next segment (to which a folded collection
     *                       corresponds).
     * @param currWeight     the weight of the count vector corresponding to the
     *                       current collection.
     * @param minWeight      the minimum weight of a solution for the entire input.
     * @param sigOracle      a signature oracle (for an informed search speedup).
     * @return A list of foldings of the current collection. The weight of the
     * count vector corresponding to each returned folding is at least
     * {@code minWeight}, and its signature lexicographic rank is small enough
     * to allow future foldings that yield solutions to the entire input.
     */
    public void allFoldings(List<PalindromeCollection> foldings, TDoubleList foldingWeights,
                            Weights w, int l, double currWeight, double minWeight, SigCurve sigOracle) {
        int minSize = w.getHeaviestCount(l);
        double currMinWeight = minWeight / currWeight;
        while (minSize > 1 && w.getWeight(l, minSize - 1) >= currMinWeight) {
            --minSize;
        }

        TDoubleList sizeWeights = new TDoubleArrayList();
        for (int i = minSize; w.getWeight(l, i) >= currMinWeight; ++i) {
            sizeWeights.add(w.getWeight(l, i) * currWeight);
        }

        List<Block> blocks = new ArrayList<>();
        for (int i = 0; i < uniqueSize(); ++i) {
            BFBPalindrome p = get(i);
            if (p instanceof Block) {
                blocks.add((Block) p);
            }
        }

        allFoldingsSelectCenters(foldings, foldingWeights, minSize, sizeWeights,
                sigOracle, blocks, blocks.size() - 1, getFirst(), 0, new Signature(), false,
                new PalindromeCollection());
    }

    private void allFoldingsSelectCenters(List<PalindromeCollection> foldings, TDoubleList foldingWeights,
                                          int minSize, TDoubleList sizeWeights, SigCurve sigOracle,
                                          List<Block> blocks, int wrapIx, BFBPalindrome center, int minCenterDeg,
                                          Signature sig, boolean centerDegIncreased, PalindromeCollection expendable) {

        if (wrapIx < 0) {
            // Base case: no additional composite palindromes may be generated.
            // Generating all foldings obtained by adding empty palindromes to
            // this collection:
            int emptyCount = 0;
            if (fullSize < minSize) {
                emptyCount = minSize - fullSize;
                add(EmptyPalindrome.SINGLETON, emptyCount);
            }

            int maxSize = minSize + sizeWeights.size() - 1;
            for (; fullSize <= maxSize; ++emptyCount) {
                sig.accommodate(this);
                if (sigOracle.withinValidRange(sig, sizeWeights.get(fullSize - minSize))) {
                    //					if (!foldings.contains(this)){
                    foldings.add(new PalindromeCollection(this));
                    foldingWeights.add(sizeWeights.get(fullSize - minSize));
                    //					}
                }
                add(EmptyPalindrome.SINGLETON, 1);
            }
            removeUnique(EmptyPalindrome.SINGLETON, emptyCount);
        } else {
            // Adding recursively additional foldings, which include concatenated
            // elements from this collection.

            Block wrap = blocks.get(wrapIx);
            int wrapMultiplicity = getMultiplicity(wrap);

            // Advancing center to the next relevant element, if needed.
            for (; center != null && (center.depth() > wrap.depth()
                    || (center.depth() == wrap.depth() && center.centerDeg() !=
                    minCenterDeg)); center = getNext(center))
                ;

            if (center != null && minCenterDeg > 0 && center.depth() < wrap.depth()) {
                center = null;
            }

            if (center == null) {
                // No more nonempty centers can participate in composite
                // palindromes.

                allFoldingsExpendComposite(foldings, foldingWeights, minSize,
                        sizeWeights, sigOracle, blocks, wrapIx, minCenterDeg,
                        sig, centerDegIncreased, expendable, getFirst(),
                        (CompositPalindrome) expendable.getFirst());

                // If minCenterDeg == 0, empty centers can be added.
                if (minCenterDeg == 0) {
                    int maxCenters = wrapMultiplicity / 2;
                    CompositPalindrome emptyCenterComposite = CompositPalindrome.make(
                            ConvexPalindrome.make(EmptyPalindrome.SINGLETON, null), wrap);
                    for (int i = 0; i < maxCenters; ++i) {
                        expendable.add(emptyCenterComposite);
                        removeUnique(wrap, 2);
                        allFoldingsExpendComposite(foldings, foldingWeights, minSize,
                                sizeWeights, sigOracle, blocks, wrapIx, minCenterDeg,
                                sig, true, expendable, getFirst(),
                                (CompositPalindrome) expendable.getFirst());
                    }
                    add(wrap, maxCenters * 2);
                    expendable.removeUnique(emptyCenterComposite, maxCenters);
                }
            } else {
                // Solving recursively without using the current center:
                allFoldingsSelectCenters(foldings, foldingWeights, minSize,
                        sizeWeights, sigOracle, blocks, wrapIx, getNext(center),
                        minCenterDeg, sig, centerDegIncreased, expendable);

                // Solving recursively with the current center used at least once:
                int maxCenters;
                if (center.equals(wrap)) {
                    // Each palindrome contains three occurrences of center/wrap.
                    maxCenters = wrapMultiplicity / 3;
                } else {
                    // Each palindrome contains two occurrences of wrap and one of center.
                    maxCenters = Math.min(getMultiplicity(center), wrapMultiplicity / 2);
                }
                CompositPalindrome composite = CompositPalindrome.make(
                        ConvexPalindrome.make(center, null), wrap);
                for (int i = 0; i < maxCenters; ++i) {
                    expendable.add(composite);
                    removeUnique(wrap, 2);
                    removeUnique(center, 1);

                    allFoldingsSelectCenters(foldings, foldingWeights, minSize,
                            sizeWeights, sigOracle, blocks, wrapIx,
                            getNext(center), minCenterDeg, sig, true, expendable);

                }
                expendable.removeUnique(composite, maxCenters);
                add(wrap, maxCenters * 2);
                add(center, maxCenters);
            }
        }
    }

    private void allFoldingsExpendComposite(List<PalindromeCollection> foldings, TDoubleList foldingWeights,
                                            int minSize, TDoubleList sizeWeights, SigCurve sigOracle,
                                            List<Block> blocks, int wrapIx, int minCenterDeg,
                                            Signature sig, boolean centerDegIncreased,
                                            PalindromeCollection expandable, BFBPalindrome expending, CompositPalindrome toExpand) {

        if (toExpand == null) {
            // No further foldings may be done using the current wrap.
            add(expandable);

            // Validity check: if the lexicographic rank of the current signature
            // is already too high, this collection can be discarded (any
            // consecutive foldings may only increase the signature).
            sig.accommodate(this);
            if (sigOracle.withinValidRange(sig)) {
                int wrapDepth = blocks.get(wrapIx).depth();
                if (wrapIx == 0 || blocks.get(wrapIx - 1).depth() > wrapDepth) {
                    // The current wrap is first among blocks of its depth, and
                    // it might be possible to generate palindromes with higher
                    // center degree.
                    if (centerDegIncreased) {
                        // Finding the last block of the same depth:
                        int lastInDepthIx = wrapIx + 1;
                        for (; lastInDepthIx < blocks.size() && blocks.get(lastInDepthIx).depth()
                                == wrapDepth; ++lastInDepthIx)
                            ;
                        --lastInDepthIx;
                        allFoldingsSelectCenters(foldings, foldingWeights, minSize, sizeWeights, sigOracle,
                                blocks, lastInDepthIx, getFirst(), minCenterDeg + 1, sig, false, new PalindromeCollection());
                    } else {
                        allFoldingsSelectCenters(foldings, foldingWeights, minSize, sizeWeights, sigOracle,
                                blocks, wrapIx - 1, getFirst(), 0, sig, false, new PalindromeCollection());
                    }
                } else {
                    allFoldingsSelectCenters(foldings, foldingWeights, minSize, sizeWeights, sigOracle,
                            blocks, wrapIx - 1, getFirst(), minCenterDeg, sig, centerDegIncreased, new PalindromeCollection());
                }
            }
            remove(expandable);
        } else {
            for (; expending != null && expending.depth() >= toExpand.center.minPartDepth();
                 expending = getNext(expending))
                ;
            if (!expandable.contains(toExpand) || expending == null) {
                // All occurrences of toExpand where already expanded and/or all
                // expansions to toExpand were done..
                allFoldingsExpendComposite(foldings, foldingWeights, minSize,
                        sizeWeights, sigOracle, blocks, wrapIx, minCenterDeg,
                        sig, centerDegIncreased, expandable, getFirst(),
                        (CompositPalindrome) expandable.getNext(toExpand));
            } else if (!contains(expending)) {
                throw new RuntimeException("Expansion bug!!!");
            } else {
                // Solving recursively without using the current expending element:
                allFoldingsExpendComposite(foldings, foldingWeights, minSize, sizeWeights,
                        sigOracle, blocks, wrapIx, minCenterDeg, sig,
                        centerDegIncreased, expandable, getNext(expending), toExpand);

                // Expanding 'toExpand' using 'expanding':
                int nestingDeg = ((ConvexPalindrome) toExpand.center).nestingDeg();
                int toExpandMultiplicity = expandable.getMultiplicity(toExpand);
                int expandingMultiplicity = getMultiplicity(expending);
                int expansionFactor = 1 << nestingDeg;
                int maxExpentions = Math.min(expandingMultiplicity / expansionFactor,
                        toExpandMultiplicity);
                CompositPalindrome composite = ((CompositPalindrome) toExpand).extendCenter(expending);

                for (int i = 0; i < maxExpentions; ++i) {
                    expandable.add(composite);
                    expandable.removeUnique(toExpand, 1);
                    removeUnique(expending, expansionFactor);

                    allFoldingsExpendComposite(foldings, foldingWeights, minSize,
                            sizeWeights, sigOracle, blocks, wrapIx, minCenterDeg,
                            sig, true, expandable,
                            getNext(expending), toExpand);
                }

                add(expending, maxExpentions * expansionFactor);
                expandable.add(toExpand, maxExpentions);
                expandable.removeUnique(composite, maxExpentions);
            }
        }
    }

    public boolean allFoldings1(FoldingHandler foldingHandler, double currWeight,
                                Weights w, int l, double minWeight, SigCurve sigOracle, Signature sig) {

        return allFoldingsSelectCenters1(foldingHandler, w, l, currWeight, minWeight, sigOracle, lastBlock(),
                getFirst(), 0, sig, false, new PalindromeCollection());
    }

    private Block lastBlock() {
        for (int i = uniqueSize() - 1; i >= 0; --i) {
            BFBPalindrome p = get(i);
            if (p instanceof Block) {
                return (Block) p;
            }
        }
        return null;
    }


    //	-------------------------------------------------------------------

    private boolean allFoldingsSelectCenters1(FoldingHandler foldingHandler,
                                              Weights w, int l, double currWeight, double minWeight,
                                              SigCurve sigOracle, Block currWrap, BFBPalindrome currCenter,
                                              int currCenterDeg, Signature sig, boolean centerDegIncreased,
                                              PalindromeCollection expandable) {

        boolean continueEnumeration = true;
        if (currWrap == null) {
            // Base case: no additional composite palindromes may be generated.
            // Generating all foldings obtained by adding empty palindromes to
            // this collection:
            int emptyCount = 0;
            double nextMinWeight = minWeight / currWeight;
            int minSize = w.getMinCount(l, nextMinWeight);
            int maxSize = w.getMaxCount(l, nextMinWeight);

            if (fullSize < minSize) {
                emptyCount = minSize - fullSize;
                add(EmptyPalindrome.SINGLETON, emptyCount);
            }

            for (; continueEnumeration && fullSize <= maxSize; ++emptyCount) {
                sig.accommodate(this);
                double weight = currWeight * w.getWeight(l, fullSize);
                if (sigOracle.withinValidRange(sig, weight)) {
                    continueEnumeration = foldingHandler.handle(this, l - 1, weight);
                }
                add(EmptyPalindrome.SINGLETON, 1);
            }
            removeUnique(EmptyPalindrome.SINGLETON, emptyCount);
        } else {
            // Adding recursively additional foldings, which include concatenated
            // elements from this collection.

            int wrapMultiplicity = getMultiplicity(currWrap);

            // Advancing center to the next relevant element, if needed.
            for (; currCenter != null && (currCenter.depth() > currWrap.depth()
                    || (currCenter.depth() == currWrap.depth() && currCenter.centerDeg() !=
                    currCenterDeg)); currCenter = getNext(currCenter))
                ;

            if (currCenter != null && currCenterDeg > 0 && currCenter.depth() < currWrap.depth()) {
                currCenter = null;
            }

            if (currCenter == null) {
                // No more nonempty centers can participate in composite
                // palindromes.

                continueEnumeration = allFoldingsExpendComposite1(foldingHandler, w, l, currWeight, minWeight,
                        sigOracle, currWrap, currCenterDeg,
                        sig, centerDegIncreased, expandable, getFirst(),
                        (CompositPalindrome) expandable.getFirst());

                // If currCenterDeg == 0, empty centers can be added.
                if (currCenterDeg == 0) {
                    int maxCenters = wrapMultiplicity / 2;
                    CompositPalindrome emptyCenterComposite = CompositPalindrome.make(
                            ConvexPalindrome.make(EmptyPalindrome.SINGLETON, null), currWrap);
                    int i = 0;
                    for (; continueEnumeration && i < maxCenters; ++i) {
                        expandable.add(emptyCenterComposite);
                        removeUnique(currWrap, 2);
                        continueEnumeration = allFoldingsExpendComposite1(foldingHandler, w, l, currWeight,
                                minWeight, sigOracle, currWrap, currCenterDeg,
                                sig, true, expandable, getFirst(),
                                (CompositPalindrome) expandable.getFirst());
                    }
                    add(currWrap, i * 2);
                    expandable.removeUnique(emptyCenterComposite, i);
                }
            } else {
                // Solving recursively without using the current center:
                continueEnumeration = allFoldingsSelectCenters1(foldingHandler, w, l, currWeight,
                        minWeight, sigOracle, currWrap, getNext(currCenter),
                        currCenterDeg, sig, centerDegIncreased, expandable);

                // Solving recursively with the current center used at least once:
                int maxCenters;
                if (currCenter.equals(currWrap)) {
                    // Each palindrome contains three occurrences of center/wrap.
                    maxCenters = wrapMultiplicity / 3;
                } else {
                    // Each palindrome contains two occurrences of wrap and one of center.
                    maxCenters = Math.min(getMultiplicity(currCenter), wrapMultiplicity / 2);
                }

                CompositPalindrome composite = CompositPalindrome.make(
                        ConvexPalindrome.make(currCenter, null), currWrap);
                int i = 0;
                for (; continueEnumeration && i < maxCenters; ++i) {
                    expandable.add(composite);
                    removeUnique(currWrap, 2);
                    removeUnique(currCenter, 1);

                    continueEnumeration = allFoldingsSelectCenters1(foldingHandler, w, l, currWeight,
                            minWeight, sigOracle, currWrap, getNext(currCenter),
                            currCenterDeg, sig, true, expandable);
                }
                expandable.removeUnique(composite, i);
                add(currWrap, i * 2);
                add(currCenter, i);
            }
        }
        return continueEnumeration;
    }

    private boolean allFoldingsExpendComposite1(FoldingHandler foldingHandler,
                                                Weights w, int l, double currWeight, double minWeight,
                                                SigCurve sigOracle, Block currWrap, int currCenterDeg,
                                                Signature sig, boolean centerDegIncreased, PalindromeCollection expandable,
                                                BFBPalindrome expending, CompositPalindrome toExpand) {

        boolean continueEnumeration = true;
        if (toExpand == null) {
            // No further foldings may be done using the current wrap.
            add(expandable);

            // Validity check: if the lexicographic rank of the current signature
            // is already too high, this collection can be discarded (any
            // consecutive foldings may only increase the signature).
            sig.accommodate(this);
            if (sigOracle.withinValidRange(sig)) {
                int wrapDepth = currWrap.depth();
                Block prevWrap = prevBlock(currWrap);
                if (prevWrap == null || prevWrap.depth() > wrapDepth) {
                    // The current wrap is first among blocks of its depth, and
                    // it might be possible to generate palindromes with higher
                    // center degree.
                    if (centerDegIncreased) {
                        // Finding the last block of the same depth:
                        Block lastInDepth = lastInDepth(currWrap);
                        if (lastInDepth != null) {
                            continueEnumeration = allFoldingsSelectCenters1(foldingHandler, w, l, currWeight,
                                    minWeight, sigOracle, lastInDepth, getFirst(),
                                    currCenterDeg + 1, sig, false, new PalindromeCollection());
                        } else {
                            continueEnumeration = allFoldingsSelectCenters1(foldingHandler, w, l, currWeight,
                                    minWeight, sigOracle, prevWrap, getFirst(),
                                    0, sig, false, new PalindromeCollection());
                        }
                    } else {
                        continueEnumeration = allFoldingsSelectCenters1(foldingHandler, w, l, currWeight,
                                minWeight, sigOracle, prevWrap, getFirst(),
                                0, sig, false, new PalindromeCollection());
                    }
                } else {
                    continueEnumeration = allFoldingsSelectCenters1(foldingHandler, w, l, currWeight,
                            minWeight, sigOracle, prevWrap, getFirst(),
                            currCenterDeg, sig, centerDegIncreased, new PalindromeCollection());
                }
            }
            remove(expandable);
        } else {
            for (; expending != null && expending.depth() >= toExpand.center.minPartDepth();
                 expending = getNext(expending))
                ;
            if (!expandable.contains(toExpand) || expending == null) {
                // All occurrences of toExpand where already expanded and/or all
                // expansions to toExpand were done.
                continueEnumeration = allFoldingsExpendComposite1(foldingHandler, w, l, currWeight,
                        minWeight, sigOracle, currWrap, currCenterDeg,
                        sig, centerDegIncreased, expandable, getFirst(),
                        (CompositPalindrome) expandable.getNext(toExpand));
            } else if (!contains(expending)) {
                throw new RuntimeException("Expansion bug!!!");
            } else {
                // Solving recursively without using the current expending element:
                continueEnumeration = allFoldingsExpendComposite1(foldingHandler, w, l, currWeight,
                        minWeight, sigOracle, currWrap, currCenterDeg,
                        sig, centerDegIncreased, expandable, getNext(expending),
                        toExpand);

                // Expanding 'toExpand' using 'expanding':
                int nestingDeg = ((ConvexPalindrome) toExpand.center).nestingDeg();
                int toExpandMultiplicity = expandable.getMultiplicity(toExpand);
                int expandingMultiplicity = getMultiplicity(expending);
                int expansionFactor = 1 << nestingDeg;
                int maxExpentions = Math.min(expandingMultiplicity / expansionFactor,
                        toExpandMultiplicity);
                CompositPalindrome composite = ((CompositPalindrome) toExpand).extendCenter(expending);

                int i = 0;
                for (; continueEnumeration && i < maxExpentions; ++i) {
                    expandable.add(composite);
                    expandable.removeUnique(toExpand, 1);
                    removeUnique(expending, expansionFactor);

                    continueEnumeration = allFoldingsExpendComposite1(foldingHandler, w, l, currWeight,
                            minWeight, sigOracle, currWrap, currCenterDeg,
                            sig, true, expandable, getNext(expending),
                            toExpand);
                }

                add(expending, i * expansionFactor);
                expandable.add(toExpand, i);
                expandable.removeUnique(composite, i);
            }
        }
        return continueEnumeration;
    }

    private Block lastInDepth(Block b) {
        Block lastInDepth = null;
        int ix = Collections.binarySearch(this, b);
        if (ix < 0) {
            ix = -ix - 1;
        } else {
            lastInDepth = b;
            ++ix;
        }
        int size = uniqueSize();
        int depth = b.depth();
        for (; ix < size && get(ix).depth() == depth; ++ix) {
            if (get(ix) instanceof Block) {
                lastInDepth = (Block) get(ix);
            }
        }
        return lastInDepth;
    }

    private Block prevBlock(Block block) {
        int ix = Collections.binarySearch(this, block);
        if (ix < 0) {
            ix = -ix - 1;
        }
        --ix;
        for (; ix >= 0 && !(get(ix) instanceof Block); --ix) ;
        Block b = null;
        if (ix >= 0) b = (Block) get(ix);
        return b;
    }

    private BFBPalindrome getNext(BFBPalindrome element) {
        int ix = Math.abs(Collections.binarySearch(this, element) + 1);
        if (ix == uniqueSize()) return null;
        else return get(ix);
    }

    private BFBPalindrome getFirst() {
        if (isEmpty()) return null;
        else return get(0);
    }

    public int getMultiplicity(BFBPalindrome p) {
        int ix = Collections.binarySearch(this, p);
        if (ix >= 0) return getMultiplicity(ix);
        else return 0;
    }

    private void add(PalindromeCollection toAdd) {
        for (int i = toAdd.uniqueSize() - 1; i >= 0; --i) {
            add(toAdd.get(i), toAdd.getMultiplicity(i));
        }
    }

    private void remove(PalindromeCollection toRemove) {
        for (int i = toRemove.uniqueSize() - 1; i >= 0; --i) {
            removeUnique(toRemove.get(i), toRemove.getMultiplicity(i));
        }
    }

    public int getSignature(int[] s) {
        Arrays.fill(s, 0);
        if (fullSize > 0) {
            int unique = uniqueSize();
            int currMultiplicity = multiplicities.get(0);
            int d = Env.dig(currMultiplicity);
            s[d] = 1;
            s[d + 1] = -currMultiplicity / (1 << (d + 1));
            int r = d + 1;

            for (int i = 1; i < unique; ++i) {
                currMultiplicity = multiplicities.get(i);
                d = Env.dig(currMultiplicity);
                if (d >= r && get(i).depth() == get(i - 1).depth()) {
                    s[r] -= currMultiplicity / (1 << r);
                } else {
                    int sum = 0;
                    for (int j = r; j > d; --j) {
                        sum = sum * 2 + Math.abs(s[j]);
                        s[j] = 0;
                    }
                    if (s[d] < 0) {
                        ++sum;
                    }
                    ++s[d];
                    s[d + 1] = -sum - currMultiplicity / (1 << (d + 1));
                    r = d + 1;
                }
            }
            return r;
        } else return 0;
    }

    public BFBPalindrome removeUnique(int ix, int multiplicity) {
        if (multiplicity > 0) {
            if (ix < 0 || ix >= uniqueSize()) {
                throw new IllegalArgumentException("Illegal unique element index: "
                        + ix + " (the collection contains " + uniqueSize() +
                        " unique elements).");
            }

            int currMultiplicity = multiplicities.get(ix);
            BFBPalindrome p;
            if (currMultiplicity > multiplicity) {
                multiplicities.set(ix, currMultiplicity - multiplicity);
                p = get(ix);
                fullSize -= multiplicity;
            } else if (currMultiplicity == multiplicity) {
                p = super.remove(ix);
                multiplicities.remove(ix);
                fullSize -= currMultiplicity;
            } else {
                throw new IllegalArgumentException("Number of copies to remove (" +
                        multiplicity + ") is greater than the number of copies " +
                        "in the collection (" + currMultiplicity + ").");
            }
            return p;
        } else return null;
    }

    public Palindrome removeUnique(BFBPalindrome p, int multiplicity) {
        return removeUnique(Collections.binarySearch(this, p), multiplicity);
    }

    public int fullSize() {
        return fullSize;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = super.hashCode();
        result = prime * result + fullSize;
        result = prime * result
                + ((multiplicities == null) ? 0 : multiplicities.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (!super.equals(obj))
            return false;
        if (getClass() != obj.getClass())
            return false;
        PalindromeCollection other = (PalindromeCollection) obj;
        if (fullSize != other.fullSize)
            return false;
        if (multiplicities == null) {
            if (other.multiplicities != null)
                return false;
        } else if (!multiplicities.equals(other.multiplicities))
            return false;
        return true;
    }

}

