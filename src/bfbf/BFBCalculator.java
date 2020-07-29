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


import bfbf.weights.Weights;

import java.util.*;
import java.util.Map.Entry;

import static java.lang.StrictMath.min;

/**
 * This class implements the algorithms presented in the paper "Algorithms for
 * Breakage-Fusion-Bridge detection in tumor genomes" by Shay Zakov, Marcus
 * Kinsella, and Vineet Bafna (PNAS 2013).
 *
 * @author Shay Zakov
 */
public class BFBCalculator {


    /**
     * An algorithm for the search variant of the BFB count vector problem.
     *
     * @param counts the input counts
     * @return a BFB string admitting the input counts if such a string exists,
     * and otherwise the string "FAILD".
     */
    public static String searchBFB(int[] counts) {
        if (!decisionBFB(counts)) {
            return "FAILED";
        }

        int k = counts.length;
        int maxR = Env.getMaxR(counts);
        Decomposition decomposition = new Decomposition(maxR, k + 1);

        for (int l = k - 1; l >= 0; --l) {
            decomposition.fold(counts[l]);
            decomposition.wrap();
        }

        decomposition.fold(1);
        String palindrome = decomposition.L[0].get(0).toString();
        return palindrome.substring(0, palindrome.length() / 2);
    }

    /**
     * An algorithm for the decision variant of the BFB count vector problem.
     * This algorithm does not maintains explicit l-block collections, but
     * only their implied signatures.
     *
     * @param counts the input counts
     * @return {@code true} if and only if there is a BFB string admitting
     * the input count vector.
     */
    public static boolean decisionBFB(int[] counts) {
        int maxR = Env.getMaxR(counts);
        int[] s = Env.borrowIntArray(maxR + 1); // the signature array

        int r = 0;
        int k = counts.length;
        int prevCount = 0;
        boolean res = decisionBFB(counts, s, r, 0, k, prevCount);

        Env.returnIntArray(s);
        return res;
    }

    public static boolean decisionBFB(int[] counts, int[] s, int r, int from,
                                      int to, int prevCount) {
        boolean res = true;

        for (int l = to - 1; l >= from; prevCount = counts[l], --l) {
            r = decisionFold(s, prevCount, r, counts[l]);
            if (r < 0) {
                // The folding procedure failed.
                res = false;
                break;
            }
        }

        if (res) {
            res = hasPalindromeConcatenation(s, counts[from], r); // Checking if the final (implicit) collection can be folded into a single palindrome.
        }
        return res;
    }

    /**
     * Applies an implicit collection folding, given the collection size and
     * signature.
     * <p>
     * <p>
     * Let {@code B} be an l-block collection, {@code s = s(B)} its signature,
     * and {@code n} an integer. If {@code B} can be folded into a collection
     * of size {@code n}, the procedure computes the minimum signature of such
     * a collection {@code B'} (see paper for definitions).
     *
     * @param s         The signature {@code s(B)} of the (implicit) input collection
     *                  {@code B}. This array will be modified to hold the signature {@code
     *                  s(B')} of the (implicit) output folded collection {@code B'}.
     * @param prevCount The size of {@code B}.
     * @param r         The value {@code r(B)}.
     * @param n         The size of the folded collection {@code B'}.
     * @return the value {@code r(B')} if {@code B'} exists, and otherwise
     * {@code -1}.
     */
    public static final int decisionFold(int[] s, int prevCount, int r, int n) {
        if (prevCount != n) {
            int dig = -1;
            int Delta_d = prevCount;
            int factor = 0;
            int d = -1;
            dig = Env.dig(prevCount - n); // dig = d_{prevCount - n}
            if (r < dig) {
                factor = 1 << dig; // factor = 2^dig
                // Asserting that all relevant values in the signature after
                // position r are set to 0:
                for (d = r + 1; d <= dig; ++d) {
                    s[d] = 0;
                }
                --d;
            } else {
                for (d = r + 1, factor = 1 << (r + 1); d > dig; --d, factor >>= 1) {
                    Delta_d -= (factor / 2) * Math.abs(s[d - 1]);
                }
            }

            //		Conditions sustained at this point:
            //		1. d = dig = dig(prevCount - n),
            //		2. factor = 2^d,
            //		3. Delta_d is as defined in the paper.

            if (prevCount < n) {
                ++s[d]; //Increasing the signature at position d.
                // Computing s[d+1] as implied from Claim 10, the fact that
                // r(B') = d+1 (and so B'_{d+1} and L'_{d+1} are empty), and
                // that Delta'_{d+1} =  Delta_d + factor*Math.abs(s[d]):
                s[d + 1] = (Delta_d + factor * Math.abs(s[d]) - n) / (2 * factor);
                return d + 1;
            } else {
                // Finding the maximum value for d sustaining the condition in
                // line 5 of Algorithm SIGNATURE-FOLD in the paper (SI):
                for (; d > 0 && n < factor * Math.max(s[d] + 1, 0) + Delta_d; --d, factor >>= 1) {
                    Delta_d -= (factor >> 1) * Math.abs(s[d - 1]);
                }
                if (d > 0 || n >= factor * Math.max(s[d] + 1, 0) + Delta_d) {
                    // The input collection can be folded into a collection of size n.

                    // Conditions:
                    // 1. 0 <= d <= dig is the maximum integer such that n >= factor * Math.max(s[d]+1, 0) + sumS,
                    // 2. factor = 2^d,
                    // 3. sumS = Delta_d
                    if (d == dig) {
                        ++s[d];
                    } else {
                        s[d] += 2;
                    }

                    int deltaDif = factor * Math.abs(s[d]); // Delta_{d+1} = sumS + delta_dif
                    // At this point, s is the signature of the collection B + (2^d)*epsilon or B + 2*(2^d)*epsilon.

                    if (n >= Delta_d + deltaDif) {
                        // No need to increase s[d], updating s[d+1]:
                        s[d + 1] = (Delta_d + deltaDif - n) / (2 * factor);
                    } else {
                        // There is a need to increase s[d], r <= d:
                        s[d] = (Delta_d - n) / factor;
                        s[d + 1] = 0;
                    }

                    for (r = d + 1; s[r] == 0 && s[r - 1] <= 0; --r) ;
                    return r;

                } else {
                    // The input collection cannot be folded into a collection of size n.
                    return -1;
                }
            }
        } else {
            // prevCount == n and no change is made to the signature or r.
            return r;
        }
    }

    /**
     * Decides if the input signature corresponds to an l-BFB palindrome
     * collection which may be concatenated to yield a single palindrome.
     *
     * @param s the collection's signature.
     * @param n the collection's size.
     * @param r the index of the last nonempty collection in the decomposition
     *          of the input.
     * @return {@code true} if an only if the (implicit) input collection may
     * be concatenated into a single palindrome.
     */
    public static boolean hasPalindromeConcatenation(int[] s, int n, int r) {
        if (n == 1 || s[0] == 0) {
            return true;
        } else if (s[0] > 1) {
            return false;
        } else {
            int i = 1;
            for (; i <= r && s[i] == 0; ++i) ;
            return s[i] <= 0;
        }
    }

    /**
     * An algorithm for the distance variant of the BFB count vector problem,
     * considering all substrings of the input.
     *
     * @param counts    the input count vector.
     * @param minWeight a bound over the minimum allowed approximation weight.
     * @return a matrix of doubles {@code distances}, where for i<=j
     * {@code distances[i][j]} is the distance of the sub-vector
     * {@code count[i...j]} from the closest BFB count vector.
     */
    public static double[][] distanceBFB(int[] counts, double[][] countWeights, double minWeight) {
        int k = counts.length;

        double[][] distances = new double[k][k];

        // Computing all sub-vector distances from BFB vectors:
        for (int i = counts.length; i > 0; --i) {
            Env.returnSolution(longestBFBSuffix(counts, countWeights, i, minWeight, distances));
        }

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                if (distances[i][j] < minWeight) {
                    distances[i][j] = 0;
                } else if (distances[i][j] > 1) {
                    distances[i][j] = 1; // fixing floating point artifacts.
                }
            }
        }

        return distances;
    }

    /**
     * @param counts       the input count vector.
     * @param countWeights
     * @param minWeight    a bound over the minimum allowed approximation weight.
     * @param mode         either {@code Env.SUFFIX} or {@code Env.SUBSTRING}.
     * @return the longest BFB vector which approximates a suffix/substring
     * (depending on the mode) of the input vector.
     */
    public static Solution1 longestBFB(int[] counts, double[][] countWeights, double minWeight, int mode) {
        Solution1 solution = null;

        switch (mode) {
            case Env.SUFFIX:
                solution = longestBFBSuffix(counts, countWeights, minWeight);
                break;
            case Env.SUBSTRING:
                solution = longestBFBSubstring(counts, countWeights, minWeight);
                break;
        }
        return solution;
    }

    /**
     * Computes the longest sub-vector of the input which is a BFB
     * count vector.
     *
     * @param counts       the input count vector.
     * @param countWeights
     * @param minWeight    a bound over the minimum allowed approximation weight.
     * @return a {@code Solution1} object representing a longest BFB vector
     * whose distance from some sub-vector of the input vector is at most
     * {@code minWeight}.
     */
    public static Solution1 longestBFBSubstring(int[] counts, double[][] countWeights, double minWeight) {
        Solution1 bestSolution = longestBFBSuffix(counts, countWeights, counts.length, minWeight, null);
        for (int k = counts.length - 1; k >= bestSolution.getLength(); --k) {
            Solution1 currSolution = longestBFBSuffix(counts, countWeights, k, minWeight, null);
            if ((currSolution.getLength() > bestSolution.getLength()) ||
                    ((currSolution.getLength() == bestSolution.getLength()) &&
                            (currSolution.getWeight() > bestSolution.getWeight()))) {
                Env.returnSolution(bestSolution);
                bestSolution = currSolution;
            } else {
                Env.returnSolution(currSolution);
            }
        }
        return bestSolution;
    }

    /**
     * Computes all sub-vectors of the input which are BFB
     * count vector, given a minimum sub-vector length and a
     * maximum count error.
     *
     * @param counts    the input count vector.
     * @param minWeight a bound over the minimum allowed approximation weight.
     */
    @Deprecated
    public static List<int[]> allMaximalBFBSubstring(int[] counts, int minLength, double minWeight) {
        List<int[]> solutions = new ArrayList<int[]>();

        int[] tempNeighbor = new int[counts.length];
        int[] tempSignature = new int[Env.getMaxR(counts) + 1];
        for (int to = counts.length; to >= minLength; --to) {
            Arrays.fill(tempSignature, 0);
            allBFBSubstring(counts, minLength, to - 1, to,
                    1, minWeight, solutions, tempNeighbor, tempSignature, 0, 0);
        }

        removeSubstrings(solutions, minLength);
        return solutions;
    }

    /**
     * Computes all eta-optimal sub-vectors of the input which are BFB
     * count vectors, given a minimum sub-vector length and 0 < eta <= 1.
     *
     * @param weights   count weights.
     * @param minLength the minimum sub-vector length.
     * @param eta       the minimum weight ratio between the weight of a reported
     *                  solution and the maximum weight of a count vector over the same segments.
     */

    public static List<int[]> allMaximalBFBSubstring(Weights weights, int minLength, double eta) {
        List<int[]> solutions = new ArrayList<int[]>();

        int k = weights.length();
        int[] tempNeighbor = new int[k];
        int[] tempSignature = new int[Env.getMaxR(weights) + 1];
        for (int to = k; to >= minLength; --to) {
            Arrays.fill(tempSignature, 0);
            allBFBSubstring(weights, minLength, to - 1, to,
                    1, eta, solutions, tempNeighbor, tempSignature, 0, 0);
        }

        removeSubstrings(solutions, minLength);
        return solutions;
    }

    private static void removeSubstrings(List<int[]> solutions, int minLength) {
        if (solutions.size() < 2) {
            return;
        }

        int length = solutions.get(0).length;
        ;
        final int[] first = {0};
        List<int[]> toRemove = new ArrayList<int[]>();

        Comparator<int[]> comparator = new Comparator<int[]>() {

            @Override
            public int compare(int[] o1, int[] o2) {
                int i = first[0];
                if (o1[i] == 0) {
                    return o2[i];
                } else if (o2[i] == 0) {
                    return -o1[i];
                } else {
                    for (; i < o1.length && o1[i] == o2[i]; ++i) ;
                    if (i < o1.length) {
                        return o1[i] - o2[i];
                    } else if (first[0] > 0 && o1[first[0] - 1] + o2[first[0] - 1] != 0) {
                        return o1[first[0] - 1] - o2[first[0] - 1];
                    } else return 0;
                }
            }
        };

        for (int i = 0; i < length - minLength; ++i) {
            toRemove.clear();
            first[0] = i;
            Collections.sort(solutions, comparator);
            int[] currSolution = solutions.get(0), nextSolution = solutions.get(1);
            boolean cont = nextSolution[i] > 0;
            for (int j = 2; cont; ++j) {
                if ((i == 0 || currSolution[i - 1] == 0) && isPrefix(currSolution, nextSolution, i)) {
                    toRemove.add(currSolution);
                }
                if (j < solutions.size()) {
                    currSolution = nextSolution;
                    nextSolution = solutions.get(j);
                    cont = nextSolution[i] > 0;
                } else {
                    cont = false;
                }
            }
            solutions.removeAll(toRemove);
        }
    }

    private static boolean isPrefix(int[] seq1, int[] seq2,
                                    int i) {
        for (; i < seq1.length && seq1[i] == seq2[i]; ++i) ;
        return i == seq1.length || seq1[i] == 0;
    }

    private static void allBFBSubstring(int[] counts, int minLength, int from,
                                        int to, double currError, double maxError, List<int[]> solutions,
                                        int[] tempNeighbor, int[] s, int prevCount, int r) {

        allBFBSubstring(Env.errorModel.getWeights(counts, maxError), minLength, from, to,
                currError, maxError, solutions, tempNeighbor, s, prevCount, r);
//		if (from >= 0){
//			double maxCurrError = Env.errorModel.maxCurrError(currError, maxError);
//
//			int maxFromValuue = Env.errorModel.maxRealValue(counts[from], maxCurrError);
//			int[] tmpS = Env.borrowIntArray(s.length);
//
//			for (int fromCount = Env.errorModel.minRealValue(counts[from], maxCurrError);
//					fromCount <= maxFromValuue; ++fromCount){
//				System.arraycopy(s, 0, tmpS, 0, r+1);
//				int tmpR = decisionFold(tmpS, prevCount, r, fromCount);
//				if (tmpR > 0 && hasPalindromeConcatenation(tmpS, fromCount, tmpR)){
//					tempNeighbor[from] = fromCount;
//					if (to-from >= minLength){
//						solutions.add(Arrays.copyOf(tempNeighbor, tempNeighbor.length));
//					}
//					double tmpError = Env.errorModel.accumulate(currError, fromCount, counts[from]);
//					allBFBSubstring(counts, minLength, from-1, to, tmpError, maxError, solutions, tempNeighbor, tmpS, fromCount, tmpR);
//				}
//			}
//			tempNeighbor[from] = 0;
//			Env.returnIntArray(tmpS);
//		}
    }

    private static void allBFBSubstring(Weights weights, int minLength, int from,
                                        int to, double currWeight, double minWeight, List<int[]> solutions,
                                        int[] tempNeighbor, int[] s, int prevCount, int r) {
        if (from >= 0) {
            double minCurrWeight = min(currWeight, minWeight);

            int maxFromValue = weights.getMaxCount(from, minCurrWeight);
            int[] tmpS = Env.borrowIntArray(s.length);

            for (int fromCount = weights.getMinCount(from, minCurrWeight);
                 fromCount <= maxFromValue; ++fromCount) {
                System.arraycopy(s, 0, tmpS, 0, r + 1);
                int tmpR = decisionFold(tmpS, prevCount, r, fromCount);
                if (tmpR > 0 && hasPalindromeConcatenation(tmpS, fromCount, tmpR)) {
                    tempNeighbor[from] = fromCount;
                    if (to - from >= minLength) {
                        solutions.add(Arrays.copyOf(tempNeighbor, tempNeighbor.length));
                    }
                    double tmpError = currWeight * weights.getWeight(from, fromCount);
                    allBFBSubstring(weights, minLength, from - 1, to, tmpError, minWeight,
                            solutions, tempNeighbor, tmpS, fromCount, tmpR);
                }
            }
            tempNeighbor[from] = 0;
            Env.returnIntArray(tmpS);
        }
    }

    public static Solution1 longestBFBSuffix(int[] counts, double[][] countWeights, double maxNormalizedError) {
        return longestBFBSuffix(counts, countWeights, counts.length, maxNormalizedError, null);
    }

    /**
     * Computes the longest suffix of {@code counts[0...k-1]} which is a BFB
     * count vector.
     *
     * @param counts       the input count vector.
     * @param countWeights
     * @param k            the last index to consider in {@code counts}.
     * @param minWeight    a bound over the minimum allowed approximation weight.
     * @param distances    if given, the procedure will set {@code distances[l][k]}
     *                     for every {@code 0 <= l < k} to the minimum distance of {@code
     *                     counts[l...k-1]} from a BFB count vector.
     * @return a {@code Solution1} object corresponding to a longest BFB vector
     * whose distance from a suffix of the prefix  {@code counts[0...k-1]} of the
     * input vector is at most {@code maxNormalizedError}.
     */
    public static Solution1 longestBFBSuffix(int[] counts, double[][] countWeights, int k,
                                             double minWeight, double[][] distances) {
        boolean saveDistances = distances != null;
        int maxR = Env.getMaxR(counts);
        int prevCount, r, tmpR;
        int[] s = Env.borrowIntArray(maxR + 1);
        int[] Delta = Env.borrowIntArray(maxR + 2);
        int bestLength = 0;

        // For the current loop index l, currLayerSolutions maps an integer m
        // to a list of candidate solutions to the sub-vector counts[l...k-1],
        // given that counts[l] is an obtained from m by introducing some error.
        // For every two candidate solutions for the same sub-vector, one has a
        // strictly lower distance from a BFB vector, while the other has a
        // strictly lower signature. prevLayerSolutions is defined the same with
        // respect to layer l+1.
        Map<Integer, List<Solution1>> currLayerSolutions = new HashMap<Integer, List<Solution1>>();
        Map<Integer, List<Solution1>> prevLayerSolutions = new HashMap<Integer, List<Solution1>>();
        Solution1 currSolution = Env.borrowSolution(counts.length, maxR);

        addSolution(currLayerSolutions, currSolution, 0);
        Solution1 bestSolution = Env.borrowSolution(currSolution.counts, currSolution.s, 0, currSolution.getWeight());

        for (int l = k - 1; l >= 0 && !currLayerSolutions.isEmpty(); --l) {
            Map<Integer, List<Solution1>> tmp = currLayerSolutions;
            currLayerSolutions = prevLayerSolutions;
            prevLayerSolutions = tmp;
            for (List<Solution1> list : currLayerSolutions.values()) {
                Env.returnSolutionList(list);
            }
            currLayerSolutions.clear();

            for (Entry<Integer, List<Solution1>> entry : prevLayerSolutions.entrySet()) {
                for (Solution1 solution : entry.getValue()) {
//					double maxCurrError = Env.errorModel.maxCurrError(solution.error, minWeight);
//					int maxCurrCount = Env.errorModel.maxRealValue(counts[l], maxCurrError);
//					int currCount = Env.errorModel.minRealValue(counts[l], maxCurrError);

                    double minCurrWeight = minWeight / solution.weight;
                    int i = 0;
                    if (countWeights[l].length == 0) {
                        int q = 1;
                    }
                    for (; countWeights[l][i] < minCurrWeight; ++i) ;

//					int maxCurrCount = Env.errorModel.maxRealValue(counts[l], minCurrWeight);
//					int currCount = Env.errorModel.minRealValue(counts[l], minCurrWeight);
//
//					
//					if (currCount > maxCurrCount){
//						// may happen when allowed error is small, and both observed count and maximum allowed count are 0.
//						break;
//					}
                    if (i == 0 && counts[l] == 0) {
                        ++i;
                    }
                    for (; i < countWeights[l].length && countWeights[l][i] >= minCurrWeight; ++i) {
//						for (; currCount <= maxCurrCount; ++currCount){
                        int currCount = counts[l] + i;
                        r = solution.getR();
                        System.arraycopy(solution.s.s, 0, s, 0, r + 1);
//						Env.fill(s, solution.s, r+1);
                        if (l + 1 == k) {
                            prevCount = 0;
                        } else {
                            prevCount = solution.counts[l + 1];
                        }
                        tmpR = decisionFold(s, prevCount, r, currCount);
                        if (tmpR >= 0 && hasPalindromeConcatenation(s, currCount, tmpR)) {
                            currSolution = Env.borrowSolution(solution.counts, s, tmpR,
                                    solution.weight * countWeights[l][i]);//Env.errorModel.accumulate(solution.error, currCount, counts[l]));
                            currSolution.counts[l] = currCount;
                            if (!addSolution(currLayerSolutions, currSolution, currCount)) {
                                Env.returnSolution(currSolution);
                            } else {
                                if (saveDistances && currSolution.weight > distances[l][k - 1]) { //Env.errorModel.compareErrors(currSolution.error, distances[l][k-1]) < 0){
                                    distances[l][k - 1] = currSolution.weight;
                                }
                                if (bestLength < k - l || currSolution.weight > bestSolution.weight) { //Env.errorModel.compareErrors(currSolution.error, bestSolution.error) < 0){
                                    Env.returnSolution(bestSolution);
                                    bestSolution = Env.borrowSolution(currSolution.counts, s, tmpR, currSolution.weight);
                                    bestLength = k - l;
                                }
                            }
                        }
                    }
                }
            }
        }

        Env.returnIntArray(s);
        Env.returnIntArray(Delta);
        return bestSolution;
    }


    private static boolean addSolution(Map<Integer, List<Solution1>> currLayerSolutions,
                                       Solution1 solution, int currCount) {
        List<Solution1> ls = currLayerSolutions.get(currCount);
        if (ls == null) {
            ls = Env.borrowSolutionList();
            currLayerSolutions.put(currCount, ls);
        }
        return addSolution(solution, ls);
    }

    private static boolean addSolution(Solution1 solution, List<Solution1> ls) {
        int listRank = Collections.binarySearch(ls, solution);
        if (listRank < 0) {
            listRank = -listRank - 1;

            if (listRank == 0 || solution.weight > ls.get(listRank - 1).weight) {
                while (listRank < ls.size() && solution.weight >= ls.get(listRank).weight) {
                    Env.returnSolution(ls.remove(listRank));
                }
                ls.add(listRank, solution);
            }
            return true;
        }
        return false;
    }


}
