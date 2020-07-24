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

package bfbf.weights;


/**
 * A Poisson-based error model. The Poisson probability of observing q when the
 * observation mean is r is defined by P(X = q | r) = r^q * e^{-r} / q!. The 
 * distance of q from r is defined by dist(q,r) = P(X = q | r) / P(X = r | r).
 * The distance between two vectors is the product of their entry-wise distances.
 *  
 * 
 * @author Shay Zakov
 *
 */
public class PoissonErrorModel extends ErrorModel {

	private static final double HALF_LOG_PI = Math.log(Math.PI) / 2.;
	private static final double[] factLans  = {0, 0, Math.log(2), Math.log(6), 
		Math.log(24), Math.log(120), Math.log(720), Math.log(5040), 
		Math.log(40320), Math.log(362880), Math.log(3628800), Math.log(39916800)};
	private static final int computedLans = factLans.length;

	@Override
	public double weight(int trueCount, int estimatedCount) {
		return Math.exp(poissonProbabilityApproximation(trueCount+0.5, estimatedCount)
				- poissonProbabilityApproximation(trueCount+0.5, trueCount));
	}

	public String toString(){
		return "Poisson errors";
	}

	/**
	 * Calculates an approximation of the log of a Poisson probability.
	 * 
	 * @param mean - lambda, the average number of occurrences.
	 * @param observed - the actual number of occurrences observed.
	 * @return ln(Poisson probability) - the natural log of the probability a 
	 * random variable, distributed according to the Poisson probability with
	 * parameter {@code lambda = mean}, gets the value {@code observed}.
	 */
	public static double poissonProbabilityApproximation (double mean, int observed) {
	        return observed * Math.log(mean) - mean - factorialApproximation(observed);
	}
	 
	/**
	 * Srinivasa Ramanujan ln(n!) factorial estimation.
	 * Good for larger values of n.
	 * 
	 * @param n a nonnegative integer (exception is thrown for negative values!).
	 * @return ln(n!)
	 */
	public static double factorialApproximation(int n) {
	        if (n < computedLans) return factLans[n];
	        double a = n * Math.log(n) - n;
	        double b = Math.log(n * (1. + 4. * n * (1. + 2. * n))) / 6.;
	        return a + b + HALF_LOG_PI;
	}

//	@Override
//	public void getCountWeights(int[] expectedCounts, int[] minCounts,
//			double[][] countProbabilities, double minCountProb) {
////		int[] expectedCounts = arm.getExpectedCounts();
//		double lnMinCountProb = Math.log(minCountProb);
//		
//		for (int i=0; i<expectedCounts.length; ++i){
//			double expectedCountProb = poissonProbabilityApproximation(expectedCounts[i], expectedCounts[i]);
//			int minCount = expectedCounts[i];
//			while (minCount > 0 && 
//					poissonProbabilityApproximation(minCount-1, expectedCounts[i])-expectedCountProb >= lnMinCountProb){
//				--minCount;
//			}
//
//			int maxCount = expectedCounts[i];
//			while (poissonProbabilityApproximation(maxCount+1, expectedCounts[i])-expectedCountProb >= lnMinCountProb){
//				++maxCount;
//			}
//
//			if (countProbabilities[i] == null || countProbabilities[i].length < maxCount-minCount+1){
//				countProbabilities[i] = new double[maxCount-minCount+1];
//			}
//			
//			minCounts[i] = minCount;
//			for (int count = minCount; count <= maxCount; ++count){
//				countProbabilities[i][count-minCount] = 
//						Math.exp(poissonProbabilityApproximation(count, expectedCounts[i])-expectedCountProb);
//			}
//		}
//	}

//	public void getCountWeights(int[] expectedCounts, int[] minCounts,
//			double[][] countProbabilities, double minCountProb) {
//		double lnMinCountProb = Math.log(minCountProb);
//		
//		for (int i=0; i<expectedCounts.length; ++i){
//			double lnMinWeight = lnMinCountProb + poissonProbabilityApproximation(expectedCounts[i], expectedCounts[i]);
//			int minCount = expectedCounts[i];
//			while (minCount > 0 && 
//					poissonProbabilityApproximation(minCount-1, expectedCounts[i]) >= lnMinWeight){
//				--minCount;
//			}
//
//			int maxCount = expectedCounts[i];
//			while (poissonProbabilityApproximation(maxCount+1, expectedCounts[i]) >= lnMinWeight){
//				++maxCount;
//			}
//
//			if (countProbabilities[i] == null || countProbabilities[i].length < maxCount-minCount+1){
//				countProbabilities[i] = new double[maxCount-minCount+1];
//			}
//			
//			minCounts[i] = minCount;
//			lnMinWeight -= lnMinCountProb;
//			for (int count = minCount; count <= maxCount; ++count){
//				countProbabilities[i][count-minCount] = 
//						Math.exp(poissonProbabilityApproximation(count, expectedCounts[i]) - lnMinWeight);
//			}
//		}
//	}
	
}
