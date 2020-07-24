package bfbf;

import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.Weights;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.*;


/**
 * A <b>Signature</b> is a list of integers {@code s = [s_0, s_1, ..., s_r]}, 
 * representing some organization characteristic of a BFB palindrome  
 * collection. The only restriction over element values is that the first 
 * nonzero element in a signature, if there is such, is positive. It is 
 * implicitly assumed that {@code s_d = 0} for every {@code d > r} (where 
 * {@code s_r} is the last nonzero element of the signature).  
 * <p>
 * The <b>cardinality</b> of a signature refers to the size of the collection it 
 * represents, and is given by {@code ||s|| = 2^0 * |s_0| + 2^1 * |s_1| + ... +
 * 2^r * |s_r|}. Thus, signatures can also be thought of as a generalization of 
 * binary representation of numbers (a binary representation is a signature whose 
 * elements {@code s_i} can only get the values 0 and 1).
 * <p>
 * Signatures can be sorted lexicographically. That is, for two signatures
 * {@code s = [s_0, s_1, ..., s_r]} and {@code s' = [s'_0, s'_1, ..., s'_r']}, 
 * {@code s < s'} if there is some index {@code d} such that {@code s_i = s'_i}
 * for every {@code 0 <= i < d}, and {@code s_d < s'_d}.
 * Signatures can be used in order to decide if a given collection can be folded
 * to have a specific size of elements. In particular if {@code s} is the signature 
 * of a collection {@code B}, it is possible to fold {@code B} into {@code m}
 * elements if and only if there exists a signature {@code s'} of cardinality 
 * {@code m} such that {@code s <= s'}. In this sense, the lexicographically 
 * lower the signature is, the less restrictions there are over foldings of the 
 * corresponding collection. Note that it is possible that {@code s < s'}, while 
 * the cardinality of {@code s} is bigger than the cardinality of {@code s'}. 
 *<p>
 * Some signature properties:<br>
 * 1. The number {@code a_i} of different signatures of cardinality {@code i} 
 * is given by {@code a_0 = 1, a_1 = 1}, {@code and a_i = a_{i-1} + a_{floor(i/2)}} 
 * for every {@code i > 1}. <br>
 * 2. The number {@code b_i} of different signatures of cardinality {@code i} 
 * which are lexicographically smaller than or equal to the signature 
 * {@code [1, 0, 0, ...]} is given by {@code b_i = a_{floor(i/2)}}. Therefore, 
 * {@code a_i = a_{i-1} + b_i}.
 *  
 * @author zakov
 *
 */
public class Signature implements Comparable<Signature>{

	// Static:
	//--------

	private static final Comparator<Solution1> solutionComparator = 
			new Comparator<Solution1>() {
		@Override
		public int compare(Solution1 arg0, Solution1 arg1) {
			return arg0.s.lexCompare(arg1.s);
		}
	};

	private static final int[] countComparatorStart = {0}; // an auxiliary pointer for countComparator 
	private static final int[] countComparatorEnd = {0}; // an auxiliary pointer for countComparator 

	private static final Comparator<int[]> countComparator = new Comparator<int[]>() {
		@Override
		public int compare(int[] first, int[] second) {
			for (int d = countComparatorStart[0]; d < countComparatorEnd[0]; ++d){
				if (first[d] != second[d]){
					return first[d] - second[d];
				}
			}
			return 0;
		}
	};

	public static final int L = 0;
	public static final int H = 1;
	private static int PRUNE_METHOD = 2; // for statistics experiments.


	// Non-static:
	//------------


	/**
	 * The list of signature elements
	 */
	protected TIntList s;

	/**
	 * Constructing an empty signature.
	 */
	public Signature(){
		s = new TIntArrayList();
		s.add(0);
	}

	/**
	 * Constructs the lexicographically minimal signature of cardinality {@code n}.
	 *  
	 * @param n the constructed signature's cardinality.
	 */
	public Signature(int n){
		this (n, true);
	}

	/**
	 * Constructs a signature of cardinality {@code n}, which is lexicographically 
	 * either minimal or maximal among all signatures of cardinality {@code n}.
	 *  
	 * @param n the constructed signature's cardinality.
	 * @param isMinSignature a flag signaling should a lexicographically minimal
	 * signature should be constructed. 
	 */
	public Signature(int n, boolean isMinSignature){
		s = new TIntArrayList();
		if (n < 0){
			throw new IllegalArgumentException("Carrdinality must be " +
					"nonnegative, given carrdinality: " + n +".");
		}
		else if (n > 0){
			if (isMinSignature){
				// Generating the minimum signature for a collection of size n: 
				while (n%2 == 0){
					s.add(0);
					n /= 2;
				}
				s.add(1);
				s.add(-n/2);
			}
			else{
				// Generating the maximum signature for a collection of size n: 
				s.add(n);
				s.add(0);
			}
		}
		else{
			s.add(0);
		}
	}

	
	/**
	 * Constructing a signature given the sizes of clusters in a collection 
	 * decomposition.
	 * 
	 * @param decomposition the sizes of L and H clusters in a collection 
	 * decomposition.
	 */
	public Signature(int[][] decomposition) {
		this();
		if (!isValidDecomposition(decomposition)){
			throw new IllegalArgumentException("Input decomposition is invalid.");
		}

		int r = decomposition[L].length;
		s.add(decomposition[L][0]);
		for (int i=1; i<r; ++i){
			s.add(decomposition[L][i] - decomposition[L][i-1] - 
					decomposition[H][i-1]/2 + Math.max(s.get(i-1), 0));
		}
		s.add(-decomposition[L][r-1] - 
				decomposition[H][r-1]/2 + Math.max(s.get(r-1), 0));
	}

	/**
	 * Constructing a signature out of an integer list.
	 * 
	 * @param s an integer list representing a signature.
	 */
	public Signature(List<Integer> s) {
		if (!isValidSignature(s)){
			throw new IllegalArgumentException("Input signature is invalid.");
		}
		this.s = new TIntArrayList(s.size());
		this.s.addAll(s);
	}

	/**
	 * A copy constructor.
	 * 
	 * @param toCopy the signature to copy.
	 */
	public Signature(Signature toCopy) {
		this.s = new TIntArrayList(toCopy.s);
	}

	/**
	 * 
	 * @return the cardinality of this signature.
	 */
	public int cardinality(){
		int n = 0;
		for (int d = r(); d >= 0; --d, n <<= 1){
			n += Math.abs(s.get(d));
		}
		return n >> 1;
	}


	public int r() {
		return s.size() - 1;
	}

	/**
	 * Computes sizes of decomposition collections that admit the signature.
	 * Note that there may be many decompositions corresponding to the same
	 * signature.

	 * @return the decomposition collection sizes.
	 */
	public int[][] getDecomposition(){
		int r = r();
		int[][] decomposition = new int[2][r];
		decomposition[L][r-1] = Math.max(s.get(r-1), 0) - s.get(r);
		int first = 0;
		for (; first <= r && s.get(first) == 0; ++first);
		for (int d = r-2; d>first; --d){
			decomposition[L][d] = decomposition[L][d+1] + Math.max(s.get(d), 0) - s.get(d+1);
		}

		decomposition[L][first] = s.get(first);
		if (r > first+1){
			decomposition[H][first] = 2*(decomposition[L][first+1] - s.get(first+1));
		}
		else{
			decomposition[H][first] = -2*s.get(first+1);
		}
		return decomposition;
	}

	/**
	 * Makes the minimum lexicographic increment to the signature so that it 
	 * will admit a given cardinality {@code n}.
	 * 
	 * @param n the cardinality after modification.
	 * @return {@code true} if the modification succeeded, and {@code false}
	 * otherwise. An unsuccessful modification does not change the original
	 * signature.
	 */
	public boolean minIncrement(int n) {
		return resize(n, 1);
	}

	/**
	 * Makes the minimum lexicographic decrement to the signature so that it 
	 * will admit a given cardinality {@code n}.
	 * 
	 * @param n the cardinality after modification.
	 * @return {@code true} if the modification succeeded, and {@code false}
	 * otherwise. An unsuccessful modification does not change the original
	 * signature.
	 */
	public boolean minDecrement(int n) {
		return resize(n, -1);
	}

	/**
	 * Modifying the signature so that it would match a given cardinality {@code n}.
	 * The modification makes the minimum lexicographic shift to the signature,
	 * which is either an increment or a decrement, depending on the given
	 * argument {@code increment}.
	 * 
	 * @param n the cardinality after modification.
	 * @param increment the direction of the lexicographic shift, 1 for 
	 * increasing the lexicographic rank, and -1 for decreasing the rank.
	 * The method is not guaranteed to work properly over other values.
	 * 
	 * @return {@code true} if the modification succeeded, and {@code false}
	 * otherwise. An unsuccessful modification does not change the original
	 * signature.
	 */
	private boolean resize(int n, int increment) {
		int cardinality = cardinality();
		if (n != cardinality){
			// The signature must change in order to have cardinality n.
			int d = firstPositionToModify(n, increment);
			if (d >= 0){
				int sd = s.get(d);
				int sigChange = increment;
				// The signature in position d is about to be shifted by either
				// +/-2 (if d < dig_{n-size}), or +/-1 (if d = dig_{n-size}). 
				// Due to this change, the size of the d-prefix of the signature   
				// may change. The only cases in which the prefix size doesn't 
				// change are when the value -1 at position d increases to 1, 
				// or the value 1 decreases to -1. 

				if (d < Env.dig(n-cardinality)){
					sigChange *= 2;
				}
				
				s.set(d, sd + sigChange);
				cardinality += (Math.abs(sd+sigChange) - Math.abs(sd)) << d;
				// Removing all entries greater than d in the signature:
				cardinality -= trimSuffix(d) << d;
				int toAdd = (n-cardinality) >> d; // (n-size) / 2^d

				if (toAdd < 0){
					// There is a need to further modify position d in the signature
					// by decreasing its absolute value. The sign of this position 
					// in the signature must be -increment:
					s.set(d, s.get(d) - toAdd * increment);
				}
				else{
					// Matching the signature size to n by setting the (d+1)-th 
					// element:
					s.add(-increment*toAdd/2);
				}
				trimSignature();
			}
			else{
				return false;
			}
		}
		return true;
	}

	/**
	 * 
	 * @param n new cardinality value.
	 * @param increment 1 if the signature can only increase and -1 if the 
	 * signature can only decrease.
	 * @return the maximum position d so that it is possible to modify the
	 * signature to have the new cardinality {@code n} without changing elements
	 * in positions 0...d-1. If no such position exists, -1 is returned.
	 */
	private int firstPositionToModify(int n, int increment) {
		int bestD = -1;
		int delta = 0;
		int maxD = Env.dig(cardinality() - n);
		for (int d = r(); d < maxD; ++d){
			s.add(0);
		}

		boolean isFirstNonZero = true;
		for (int d = 0, factor = 1; d <= maxD; ++d, factor *= 2){
			// invariant: factor = 2^d, delta = cardinality of the prefix of s before position d.
			int sd = s.get(d);
			if (n >= factor * Math.max(increment*sd + 1, 0) + delta){
				// Signature may be updated, but only as long as the results is
				// between [0, 0, ...] and [1, 0, ...]
				if (increment > 0){
					// Incrementing only if the result is at most [1, 0, ...]:
					if (d > 0 || (sd == 0 && maxD == 0)){
						bestD = d;	
					}
				}
				else if (!(isFirstNonZero && (sd == 0 || (sd == 1 && d < maxD)))){
					// Decrementing only if the result is at lease [0, 0, ...]:
					bestD = d;
				}
			}
			delta += factor * Math.abs(sd);
			if (sd > 0){
				isFirstNonZero = false;
			}
		}
		return bestD;
	}



	/**
	 * Removes all elements after a given position in the signature.
	 * 
	 * @param d the last position that will be kept.
	 * @return the partial contribution of the removed elements to cardinality,
	 * divided by {@code 2^d}.
	 */
	private int trimSuffix(int d) {
		int suffixSize = 0;
		for (int i = r(); i > d; --i, suffixSize *= 2){
			suffixSize += Math.abs(s.removeAt(i));
		}
		return suffixSize;
	}

	//TODO: Move to BFB/BFB_algo
	public static List<int[]> allHeavyBFBSubVectors(Weights weights, double minWeight, int minLength){

		List<int[]> heavyBFBSubVectors = new ArrayList<int[]>();
		List<int[]> toAdd = new ArrayList<int[]>();

		List<Solution1> currSolutions = new ArrayList<Solution1>();
		List<Solution1> prevSolutions = new ArrayList<Solution1>();
		List<Solution1> tmp;
		BitSet extended = new BitSet();

		Solution1 emptySolution = new Solution1();
		int[] subCounts = new int[weights.length()];
		SigCurve lSigCurves;

		for (int k = weights.length(); k >= minLength; --k){

			SigCurve[] sigCurves = lowerBoundSigCurves(weights, minWeight, k-minLength, k);
			if (sigCurves[minLength].size() == 0){
				continue;
			}

			currSolutions.add(emptySolution);
			Signature s = emptySolution.s;

			for (int l = k-1; l >= 0; --l){
				extended.clear();

				tmp = currSolutions;
				currSolutions = prevSolutions;
				prevSolutions = tmp;
				currSolutions.clear();
				int maxCount = weights.getMaxCount(l);
				double maxWeight = weights.getMaxWeight(l);

				if (k-l <= minLength){
					lSigCurves = sigCurves[minLength - k + l];
				}
				else{
					lSigCurves = null;
				}

				for (int m = weights.getMinCount(l); m <= maxCount; ++m){
					double currWeight = weights.getWeight(l, m) /  maxWeight;
					foldToCount(minWeight, currSolutions, prevSolutions,
							s, lSigCurves, m, currWeight, extended);
				}

				if (k-l > minLength){
					addSubCounts(heavyBFBSubVectors, toAdd, prevSolutions,
							extended, subCounts, k, l+1);
				}
			}

			extended.clear();
			addSubCounts(heavyBFBSubVectors, toAdd, currSolutions,
					extended, subCounts, k, 0);
		}
		return heavyBFBSubVectors;
	}

	/**
	 * Adding heavy BFB vectors which are not sub-vectors of 
	 * previous solutions, and which will not be contained
	 * in vectors added in future iterations.
	 * 
	 * @param heavyBFBSubVectors
	 * @param toAdd
	 * @param solutions
	 * @param extended
	 * @param subCounts
	 * @param k
	 * @param l
	 */
	private static void addSubCounts(List<int[]> heavyBFBSubVectors,
			List<int[]> toAdd, List<Solution1> solutions, BitSet extended,
			int[] subCounts, int k, int l) {

		toAdd.clear();
		countComparatorStart[0] = l;
		countComparatorEnd[0] = k;
		Collections.sort(heavyBFBSubVectors, countComparator);

		for (int i=extended.previousClearBit(solutions.size()-1); i >= 0; i = extended.previousClearBit(i-1)){
			Arrays.fill(subCounts, 0);
			int[] counts = solutions.get(i).counts;
			System.arraycopy(counts, 0, subCounts, l, counts.length);
			if (Collections.binarySearch(heavyBFBSubVectors, subCounts, countComparator) < 0){
				// The current counts are not contained in any previous counts in heavyBFBSubVectors
				toAdd.add(Arrays.copyOf(subCounts, subCounts.length));
			}
		}

		heavyBFBSubVectors.addAll(toAdd);
	}


	/**
	 * Finds the heaviest BFB count vector with respect to a given weight function.
	 * 
	 * @param weights a weight function.
	 * @param minLength the minimum number of segments in a solution.
	 * @param minWeight the minimum solution weight.
	 * @return A maximum weight BFB count vector with respect to the input 
	 * weights, whose length is at least {@code minLength} and whose weight
	 * it at least {@code minWeight}. If no such vector exists, null is returned.
	 */

	public static Solution1 heaviestBFBVector(Weights weights, int minLength, double minWeight){

		double bestWeight = 2;
		SigCurve[] bestCurves = null;
		int bestStart = weights.length();

		for (int start = weights.length() - minLength; start >=0; --start){
			SigCurve[] curves = lowerBoundSigCurves(weights, minWeight, start, start + minLength);
			if (curves[minLength].size() > 0){
				if (curves[minLength].weights.get(0) < bestWeight){
					bestWeight = curves[minLength].weights.get(0);
					bestCurves = curves;
					bestStart = start;
				}
			}
		}

		if (bestCurves != null){
			int[] counts = new int[weights.length()];
			for (int l = bestStart + minLength - 1; l >= bestStart; --l){

			}
		}

		else {
			return null;
		}


		List<Solution1> currSolutions = new ArrayList<Solution1>();
		List<Solution1> prevSolutions = new ArrayList<Solution1>();
		List<Solution1> tmp;
		Solution1 emptySolution = new Solution1();
		Solution1 tmpSolution = new Solution1();
		int[] heaviestCounts = weights.getHeaviestCounts();
		double currWeight;

		Solution1 optSolution = null;
		int optSolutionStart = -1;

		for (int k = weights.length(); k >= minLength; --k){
			currSolutions.clear();
			prevSolutions.clear();
			prevSolutions.add(emptySolution);

			for (int l = k-1; l >=0; --l){
				for (int i = prevSolutions.size()-1; i>=0; --i){
					Solution1 prevSolution = prevSolutions.get(i);
					currWeight = prevSolution.getWeight();
					for (int m = heaviestCounts[l]; currWeight >= minWeight; ++m, 
							currWeight = weights.getWeight(l, m) * prevSolution.getWeight()){
						extend(prevSolution, m, currWeight, k-l, tmpSolution, currSolutions);
					}

					currWeight = weights.getWeight(l, heaviestCounts[l]-1) * prevSolution.getWeight();
					for (int m = heaviestCounts[l]-1; currWeight >= minWeight; --m, 
							currWeight = weights.getWeight(l, m) * prevSolution.getWeight()){
						extend(prevSolution, m, currWeight, k-l, tmpSolution, currSolutions);
					}
				}

				if (k-l >= minLength && !currSolutions.isEmpty()){
					Solution1 currSolution = currSolutions.get(currSolutions.size()-1);
					prevSolutions.clear(); //TODO: pooling
					if (optSolution == null || 
							(currSolution.getWeight() > optSolution.getWeight()) ||
							((currSolution.getWeight() == optSolution.getWeight()) &&
									currSolution.getLength() > optSolution.getLength())){
						optSolution = currSolution;
						optSolutionStart = l;
						minWeight = optSolution.getWeight();
						prevSolutions.add(optSolution);
					}
				}
				else{
					tmp = currSolutions;
					currSolutions = prevSolutions;
					prevSolutions = tmp;
				}

				if (prevSolutions.isEmpty()){
					break;
				}
				currSolutions.clear(); //TODO: pooling
			}

		}

		if (optSolution != null){
			int[] counts = new int[weights.length()];
			System.arraycopy(optSolution.counts, 0, counts, optSolutionStart, optSolution.counts.length);
			optSolution.counts = counts;
		}

		return optSolution;
	}

	private static void extend(Solution1 prevSolution,
			int m, double currWeight, int length, Solution1 tmpSolution,
			List<Solution1> currSolutions) {
		tmpSolution.s.s.clear();
		tmpSolution.s.s.addAll(prevSolution.s.s);
		if (tmpSolution.s.minIncrement(m) && tmpSolution.s.hasPalindromicConcatenation()){
			int rank = Collections.binarySearch(currSolutions, tmpSolution, solutionComparator);
			if (rank < 0){
				rank = -rank - 1;
				if (rank == 0 || (rank > 0 && 
						currSolutions.get(rank-1).getWeight() < currWeight)){
					// A new solution is added to the list
					int[] counts = new int[length];
					counts[0] = m;
					System.arraycopy(prevSolution.counts, 0, counts, 1, length-1);
					currSolutions.add(rank, new Solution1(counts, new Signature(tmpSolution.s), currWeight));
				}
			}
			else {
				// A solution with identical signature already 
				// exists - the new solution is added only if
				// it is of higher weight. Since the signatures
				// are identical, the two solutions must have
				// the same count m for segment l. 
				Solution1 currSolution = currSolutions.get(rank);
				if (currSolution.getWeight() < currWeight){
					currSolution.setWeight(currWeight);
					System.arraycopy(prevSolution.counts, 0, currSolution.counts, 1, length-1);
				}
			}

			// removing solutions dominated by the added solution:
			while (currSolutions.size() > rank+1 && currSolutions.get(rank+1).getWeight() <= currWeight){
				currSolutions.remove(rank + 1);
			}
		}
	}

	private boolean hasPalindromicConcatenation() {
		if (get(0) < 2){
			if (get(0) == 1){
				int r = r();
				for (int d = 1; d<=r; ++d){
					if (get(d) != 0){
						return get(d) < 0;
					}
				}
			}
			return true;
		}
		else{
			return false;
		}
	}

	private static int foldToCount(double minWeight,
			List<Solution1> currSolutions, List<Solution1> prevSolutions,
			Signature s, SigCurve lSigCurves,
			int m, double currWeight, BitSet extended) {
		int size = prevSolutions.size();
		int totalExaminedSolutions = 0;

		for (int i = 0; i < size; ++i){
			Solution1 prevSolution = prevSolutions.get(i);
			double solutionWeight = prevSolution.getWeight() * currWeight;
			if (solutionWeight < minWeight){
				continue;
			}
			++totalExaminedSolutions;
			s.setTo(prevSolution.s);
			if (s.minIncrement(m)){
				//				s.setWeight(minWeight/solutionWeight);
				boolean toAdd; 

				if (PRUNE_METHOD == 0){
					toAdd = s.firstNonZero(0) > 0;
				}
				else if (PRUNE_METHOD == 1 || lSigCurves == null){
					toAdd = s.canFold(1);
				}
				else{
					toAdd = lSigCurves.withinValidRange(s, solutionWeight);
				}


				if (toAdd){
					//					++totalGenratedSolutions[0];
					Solution1 currSolution = new Solution1(addCount(m, prevSolution.counts), s, solutionWeight);
					currSolutions.add(currSolution);
					extended.set(i);
				}
			}
		}
		return totalExaminedSolutions;
	}


	private static int[] addCount(int m, int[] counts) {
		int[] newCounts = new int[counts.length+1];
		System.arraycopy(counts, 0, newCounts, 1, counts.length);
		newCounts[0] = m;
		return newCounts;
	}



	private boolean canFold(int n) {
		return canResize(n, 1);
	}

	/**
	 * Checks if there exists some collection of size n which after folding
	 * gets this signature.
	 * 
	 * @param n a collection size.
	 * @return true if and only if there is an equal or lexicographically 
	 * smaller signature than this signature, whose cardinality is n. 
	 */
	public boolean canUnfold(int n) {
		//		return canResize(n, -1);
		int size = cardinality();
		if (size == 0){
			return false;
		}
		else if (n == size){
			return true;
		}
		else{
			int f = 0;
			int r = r();
			for (; s.get(f) == 0; ++f);
			int d = Env.dig(size - n);
			return d >= f;
		}
	}

	private boolean canResize(int n, int increment) {
		int size = cardinality();
		if (n == size){
			return true;
		}
		boolean isFirstNonZero = true;
		int maxDig = Env.dig(size - n);
		int maxD = Math.min(r(), maxDig);
		int sd;
		for (int d = 0; d <= maxD; ++d, n = (n-Math.abs(sd))/2){
			// invariant: factor = 2^d, delta = Delta_d
			sd = s.get(d);
			if (n >= Math.max(increment*(sd+increment), 0)){
				if (!(isFirstNonZero && increment < 0 && (sd == 0 || (sd == 1 && d < maxDig)))){
					return true;
				}
			}
			if (sd > 0){
				isFirstNonZero = false;
			}
		}
		return false;
	}


	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((s == null) ? 0 : s.hashCode());
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
		Signature other = (Signature) obj;
		if (s == null) {
			if (other.s != null)
				return false;
		} else if (!s.equals(other.s))
			return false;
		return true;
	}

	public String toString(){
		return "n: " + cardinality() +", r: " + r() +", s: " + s.toString();
	}


	/**
	 * Checks if a given integer list represents a valid signature. A list of
	 * integers is a valid signature if it is nonempty, its first nonzero 
	 * element is positive, and its last element is non-positive.
	 * 
	 * @param s a list of integers.
	 * @return {@code true} if and only if the input string represents a valid
	 * signature.
	 */
	public static boolean isValidSignature(List<Integer> s){
		if (s == null || s.size() == 0) return false;
		else if (s.get(s.size()-1) > 0) return false;
		else {
			int d = 0, sd = s.get(0);
			for (; sd == 0 && d < s.size(); ++d, sd = s.get(d));
			return sd >= 0;
		}
	}

	/**
	 * Making sure that the signature ends with a non-positive number, and 
	 * trimming multiple zeros at the end.  
	 */
	private void trimSignature(){
		int d = r();
		int sd = s.get(d);
		if (sd >= 0){
			if (sd > 0){
				s.add(0);
			}
			else{
				for (; d > 0 && s.get(d-1) == 0; --d);
				if (d > 0 && s.get(d-1) < 0){
					s.removeAt(d);
				}
			}
		}
	}



	private static boolean isValidDecomposition(int[][] decomposition) {
		for (int i = decomposition[L].length-1; i>=0; --i){
			if (decomposition[L][i] < 0 || decomposition[H][i] < 0 ||
					(decomposition[L][i] == 0 && decomposition[H][i] > 0)){
				return false;
			}
		}
		return true;
	}




	/**
	 * Makes the minimum decrement to a signature without changing underlying collection size.
	 */
	public void decrese(){
		int d = r();
		int sd = s.get(d);
		if (sd == 0){
			// s(d-1) > 0, and it is possible to increase the signature at position d-1
			sd = s.get(d-1);
			s.set(d-1, sd-2);
			if (sd > 1){
				s.set(d, 1);
				s.add(0);
			}
			else{
				trimSignature();
			}
		}

		else{
			// sd < 0, and there is a need to decrease s(d-1)
			int sdm1 = s.get(d-1);
			int dif = 0;
			if (sdm1 < 1) dif = -1;
			else if (sdm1 > 1) dif = 1;
			s.set(d-1, sdm1-2);
			s.set(d, -sd+dif);
			s.add(0);
			trimSignature();
		}
	}

	/**
	 * Makes the minimum increment without changing underlying collection size.
	 */
	public void increse(){
		int d = r();
		int sd = s.get(d);
		if (sd < 0){
			// it is possible to increase the signature at position d
			s.set(d, sd+2);
			if (sd < -1){
				s.add(-1);
			}
			else{
				s.add(0);
			}
		}

		else{
			// sd = 0, s(d-1) > 0, and there is a need to increase s(d-2)
			sd = s.get(d-2);
			int dif = 0;
			if (sd < -1) dif = 1;
			else if (sd > -1) dif = -1;
			s.set(d-2, sd+2);
			s.set(d-1, -s.get(d-1)-dif);
			trimSignature();
		}
	}

	public boolean canDecrease(){
		int r = r();
		if (r == 0) return false;

		int d;
		for (d = 0; s.get(d) == 0; ++d);
		if (r > d+1 || s.get(d) > 1 || s.get(d+1) > 0){
			return true;
		}
		return false;
	}

	public boolean canIncrease(){
		int r = r();
		if (r == 0 || (r == 1 && s.get(1) == 0)) return false;
		return true;
	}

	/**
	 * Computes the lexicographic relation between two signatures.
	 * 
	 * @param other a signature compared to this signature.
	 * @return A negative integer if this signature precedes lexicographically 
	 * the other signature, a positive integer if the other signature 
	 * precedes lexicographically this signature, or zero if the two signatures
	 * are identical.  
	 */
	public int lexCompare(Signature other) {
		int r = r();
		int rOther = other.r();
		int maxR = Math.max(r, rOther);
		int d = 0;
		for (; d <= maxR && get(d) == other.get(d); ++d);
		return get(d) - other.get(d);
	}

	@Override
	public int compareTo(Signature other) {
		return lexCompare(other);
	}

	protected int firstNonZero(int d) {
		int r = r();
		for (; d <= r; ++d){
			int sd = s.get(d);
			if (sd != 0) {
				return sd;
			}
		}
		return 0;
	}


	/**
	 * This method computes the number of different possible extended signatures  
	 * at given sizes. An extended signature is defined the same as a regular  
	 * signature, i.e. a list of integers s_0, s_1, ..., s_r, and the only 
	 * difference is that it does not have to sustain the restriction over 
	 * regular signatures having the first nonzero element to be positive. The  
	 * size of an extended signature defined the same as for regular signatures: 
	 * 2^0*|s_0| + 2^1*|s_1| + ... 2^r * |s_r|. The number of different regular 
	 * signatures of a given size is half of the number of extended signatures 
	 * of this size (as for any extended signature whose first nonzero is 
	 * negative there is a corresponding signature whose first nonzero is 
	 * positive).  
	 * 
	 * 
	 * @param maxSize The maximum extended signature size.
	 * @param numOfExtendedSignatureLs An integer list that will be filled with  
	 * the result. If null is given, a new list is created and returned. Otherwise,  
	 * the given list is cleared and returned at the end of the run.
	 * @return the integer list numOfExtendedSignatureLs, where the i-th element 
	 * in this list is the number of different extended signature of size i.
	 */
	public static List<Long> numOfExtendedSignatureLs(int maxSize, 
			List<Long> numOfExtendedSignatureLs){
		if (numOfExtendedSignatureLs == null){
			numOfExtendedSignatureLs = new ArrayList<Long>();
		}
		else{
			numOfExtendedSignatureLs.clear();
		}

		long currCardinality = 1; // initialized to the number of signatures of an empty collection
		for (int i = 0; i <= maxSize+1; ++i){
			numOfExtendedSignatureLs.add(currCardinality);
			currCardinality += numOfExtendedSignatureLs.get((i+1)/2);
		}
		return numOfExtendedSignatureLs;
	}

	public static List<Long> numOfSignatureLs(int maxSize){
		return numOfExtendedSignatureLs(maxSize, null);
	}

	public static long numOfSignatures(int size, List<Long> numOfExtendedSignatureLs){
		if (size  == 0){
			return 1;
		}
		else{
			return numOfExtendedSignatureLs.get(size)/2;
		}
	}

	/**
	 * Computes the number of signatures of a given size that lexicographically
	 * precedes the unit signature 1, 0, 0, ...
	 *  
	 * @param size the size of the signatures.
	 * @param numOfExtendedSignatureLs the output list from the method
	 * numOfExtendedSignatureLs(...), whose input maxSize parameter must be at
	 * least
	 * @return
	 */
	public static long numOfSignaturesSmallerThanUnit(int size, List<Long> numOfExtendedSignatureLs){
		return numOfSignatures(size/2, numOfExtendedSignatureLs);
	}

	public static List<Integer> numOfBitsInSignatureRepLs(int maxSize, 
			List<Integer> numOfBitLs){
		if (numOfBitLs == null){
			numOfBitLs = new ArrayList<Integer>();
		}
		else{
			numOfBitLs.clear();
		}

		int currCardinality = 1; // initialized to the number of signatures of an empty collection
		for (int i = 0; i <= maxSize+1; ++i){
			numOfBitLs.add(currCardinality);
			currCardinality += numOfBitLs.get((i+1)/2);
		}
		return numOfBitLs;
	}


	public int get(int d) {
		if (d <= r()){
			return s.get(d);
		}
		else{
			return 0;
		}
	}

	public void clear() {
		s.clear();
		s.add(0);
	}

	/**
	 * Sets signature state to represent the signature of a given collection.
	 * Note that this operation does not modify the signature weight.
	 * 
	 * @param collection a collection whose signature will be represented by
	 * this signature object.
	 */
	public void accommodate(PalindromeCollection collection){
		clear();
		int last = collection.uniqueSize();
		int prevDepth = -1, currDepth;

		for (int i=0; i<last; ++i){
			currDepth = collection.get(i).depth();
			add(collection.getMultiplicity(i), currDepth != prevDepth);
			prevDepth = currDepth;
		}

	}

	/**
	 * Updating the signature due to adding a new element of minimum depth to  
	 * the collection.
	 * 
	 * @param multiplicity the multiplicity of the added element.
	 * @param isLower true if and only if the added element is of strictly
	 * lower depth than the previous minimum depth element in the collection. 
	 */
	private void add(int multiplicity, boolean isLower) {
		int d = Env.dig(multiplicity);
		int suffixSum = 0;
		int r = s.size()-1;

		if (d >= r && !isLower){
			// No increment is made to the signature, and the r-th element
			// is decreased.
			suffixSum = Math.abs(s.removeAt(r));
		}
		else{
			// The d-th signature position is increased by 1, r is set to
			// d+1, and the signature at position r is accommodated to make
			// signature and collection size identical.

			if (d > r){
				// In case d > r:
				s.fill(r+1, d+1, 0);
			}
			else{
				// In case d < r:
				for (; r > d; --r){
					suffixSum += suffixSum + Math.abs(s.removeAt(r));
				}
			}

			int sd = s.get(d);
			s.set(d, sd + 1);
			if (sd < 0){
				++suffixSum;
			}
			r = d+1;
		}

		suffixSum += multiplicity >> r; // dividing multiplicity by 2^r
			s.add(-suffixSum);
	}


	public void setTo(Signature other) {
		s.clear();
		s.addAll(other.s);
	}

	public void setTo(int[] other) {
		s.clear();
		s.addAll(other);
	}

	public static SigCurve[] lowerBoundSigCurves(Weights w, double minWeight, 
			int start, int end){

		SigCurve[] curves = new SigCurve[end-start+1];
		Signature s = new Signature();
		Signature prevSig;
		SigCurve currCurve, prevCurve = new SigCurve(minWeight);
		prevCurve.sigs.add(new Signature(1));
		prevCurve.weights.add(minWeight);
		double prevWeight, countWeight, currWeight;
		int n, pointsToRemove;
		curves[0] = prevCurve;

		for (int l=start; l<end; ++l){
			currCurve = new SigCurve(minWeight);

			for (int pointIx = 0; pointIx < prevCurve.size(); ++pointIx){
				prevWeight = prevCurve.weights.get(pointIx);
				prevSig = prevCurve.sigs.get(pointIx);
				n = w.getMinCount(l, prevWeight);
				countWeight = w.getWeight(l, n);

				for (; countWeight >= prevWeight; ++n, countWeight = w.getWeight(l, n)){
					s.setTo(prevSig);
					if (s.minDecrement(n)){
						currWeight = prevWeight / countWeight;
						int insertionIx = currCurve.weights.binarySearch(currWeight);
						if (insertionIx < 0){
							insertionIx = -insertionIx-1;
							if (insertionIx > 0 && s.compareTo(currCurve.sigs.get(insertionIx-1)) <= 0){
								// The point is dominated by another point on the curve
								continue;
							}
							else{
								currCurve.add(insertionIx, new Signature(s), currWeight);
							}
						}
						else if (s.compareTo(currCurve.sigs.get(insertionIx)) > 0){
							currCurve.sigs.get(insertionIx).setTo(s);
						}
						else{
							continue;
						}

						int i = insertionIx + 1;
						for (pointsToRemove = 0; i < currCurve.size() && s.compareTo(currCurve.sigs.get(i)) >= 0;
								++i, ++pointsToRemove);
						if (pointsToRemove > 0){
							currCurve.weights.remove(insertionIx+1, pointsToRemove);
							for(; i < currCurve.size(); ++i){
								currCurve.sigs.set(i-pointsToRemove, currCurve.sigs.get(i));
							}
							for (--i; pointsToRemove > 0; --i, --pointsToRemove){
								currCurve.sigs.remove(i);
							}
						}
					}
				}
			}

			curves[l-start+1] = currCurve;
			prevCurve = currCurve;
		}

		return curves;

	}

/*
	private static class SigCurve {
		private List<Signature> sigs;
		private TDoubleList weights;


		public SigCurve() {
			sigs = new ArrayList<Signature>();
			weights = new TDoubleArrayList();
		}

		public void add(int insertionIx, Signature signature, double weight) {
			sigs.add(insertionIx, signature);
			weights.insert(insertionIx, weight);
		}

		public int size() {
			return sigs.size();
		}
	}
*/

}
