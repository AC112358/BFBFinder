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
 * An abstract error model.
 *
 * @author Shay Zakov
 */
public abstract class ErrorModel {

    private static int MAX_RATIO = 2;

    /**
     * Computing the weight of an estimated count, given a true count.
     *
     * @param trueCount      the number of times a segment repeats in a chromosome.
     * @param estimatedCount the number of times the segment was estimated to
     *                       repeat in the chromosome.
     * @return a weight reflecting the probability to call the (possibly wrong)
     * estimated count, given that the true count is the correct one.
     */
    public abstract double weight(int trueCount, int estimatedCount);

    /**
     * Returns a Weights object for a count estimation vector.
     *
     * @param counts    an integer vector with segment count estimations.
     * @param minWeight a bound to the
     * @return
     */

    public Weights getWeights(int[] counts, double minWeight) {
        Weights w = new Weights(counts.length);
        for (int i = 0; i < counts.length; ++i) {
            int minValue = counts[i];
            int maxValue = counts[i];
            if (counts[i] >= 0) {
                minValue = minRealValue(counts[i], minWeight);
                maxValue = maxRealValue(counts[i], minWeight);
            }
            double[] countWeights = new double[maxValue - minValue + 1];
            for (int j = minValue; j <= maxValue; ++j) {
                if (counts[i] >= 0) {
                    countWeights[j - minValue] = weight(j, counts[i]);
                }else{
                    countWeights[j - minValue] = 1;
                }
            }
            w.setWeights(i, minValue, countWeights);
        }
        return w;
    }


    /**
     * Computing a lower bound over a true segment count, given an estimated
     * count and a minimum estimation weight.
     *
     * @param estimatedCount an (experimentally) estimated segment count.
     * @param minWeight      a minimum estimation weight.
     * @return the minimum true count such that the estimated count weight is
     * at least the given minimum weight.
     */
    public int minRealValue(int estimatedCount, double minWeight) {
        int deltaMax = estimatedCount, deltaMin = 0, deltaMid;

        // Invariant: estimatedCount - deltaMax is an invalid true count,
        // estimatedCount - deltaMin is a valid true count.
        while (deltaMax - deltaMin > 1) {
            deltaMid = (deltaMax + deltaMin) / 2;
            if (weight(estimatedCount - deltaMid, estimatedCount) >= minWeight) {
                deltaMin = deltaMid;
            } else deltaMax = deltaMid;
        }
        return estimatedCount - deltaMin;
    }

    /**
     * Computing an upper bound over a true segment count, given an estimated
     * count and a minimum estimation weight.
     *
     * @param estimatedCount an (experimentally) estimated segment count.
     * @param minWeight      a minimum estimation weight.
     * @return the maximum true count such that the estimated count weight is
     * at least the given minimum weight.
     */
    public int maxRealValue(int estimatedCount, double minWeight) {
        int deltaMax = 1, deltaMin = 0, deltaMid;

        // Invariant: estimatedCount + deltaMin is a valid true count.
        while (weight(estimatedCount + deltaMax, estimatedCount) >= minWeight) {
            deltaMin = deltaMax;
            deltaMax *= 2;
        }

        // Invariant: estimatedCount + deltaMax is an invalid true count,
        // estimatedCount + deltaMin is a valid true count.
        while (deltaMax - deltaMin > 1) {
            deltaMid = (deltaMax + deltaMin) / 2;
            if (weight(estimatedCount + deltaMid, estimatedCount) >= minWeight) {
                deltaMin = deltaMid;
            } else deltaMax = deltaMid;
        }

        return estimatedCount + deltaMin;
    }


    //DEPRECATED:

//	@Deprecated
//	public abstract double error(int realValue, int observedValue);
//	@Deprecated
//	public abstract double accumulate (double accumulatedError, double newError);
//	@Deprecated
//	public abstract double maxCurrError(double accumulatedError, double maxError);
//	@Deprecated
//	public double initialError(){return 0;}
//	@Deprecated
//	public int compareErrors(double error1, double error2){
//		if (error1 < error2) return -1;
//		else if (error1 > error2) return 1;
//		else return 0;
//	}
//
//	@Deprecated
//	public void getCountProbabilities(ChromosomeArm arm, int[] minCounts, 
//			double[][] countProbabilities, double minCountProb){
//		getCountProbabilities(arm.getExpectedCounts(), minCounts, countProbabilities, minCountProb);
//	}
//
//	@Deprecated
//	public void getCountProbabilities(int[] expectedCounts, int[] minCounts, 
//			double[][] countProbabilities, double minCountProb){
//		getCountWeights(expectedCounts, minCounts, countProbabilities, minCountProb);
//		for (int i=0; i<minCounts.length; ++i){
//			double maxProb = 0;
//			for (int j=0; j<countProbabilities[i].length; ++j){
//				maxProb = Math.max(maxProb, countProbabilities[i][j]);
//			}
//			for (int j=0; j<countProbabilities[i].length; ++j){
//				countProbabilities[i][j] /= maxProb;
//			}
//		}
//	}

//	@Deprecated
//	public abstract void getCountWeights(int[] counts, int[] minCounts, 
//			double[][] countProbabilities, double minCountProb);

//	/**
//	 * Normalize the error to the 0-1 range, 0 means identical vectors and 1
//	 * means extremely different vectors.
//	 * 	
//	 * @param error the unnormalized error.
//	 * @return the normalized value of the input error.
//	 */
//
//	@Deprecated
//	public abstract double normlizeError(double error);
//
//	/**
//	 * The reverse of the normalizeError method.
//	 * 
//	 * @param normlizedError
//	 * @return
//	 */
//	@Deprecated
//	public abstract double deNormlizeError(double normlizedError);

//	@Deprecated
//	public int minObservedValue(int realValue, double maxError){
//		int deltaMax=1, deltaMin=0, deltaMid;
//		while (deltaMax < realValue && compareErrors(error(realValue, realValue - deltaMax), maxError) <= 0){
//			deltaMin = deltaMax;
//			if (deltaMax*2 <= realValue){
//				deltaMax *= 2;
//			}
//			else {
//				deltaMax += (realValue-deltaMax+1)/2;
//			}
//		}
//		while(deltaMax - deltaMin > 1){
//			deltaMid = (deltaMax + deltaMin)/2;
//			if (compareErrors(error(realValue, realValue - deltaMid), maxError) <= 0){
//				deltaMin = deltaMid;
//			}
//			else deltaMax = deltaMid;
//		}
//
//		return realValue - deltaMin;
//	}
//
//	@Deprecated
//	public int maxObservedValue(int realValue, double maxError){
//		int deltaMax=1, deltaMin=0, deltaMid;
//		while (deltaMax < realValue*MAX_RATIO && compareErrors(error(realValue + deltaMax, realValue), maxError) <= 0){
//			deltaMin = deltaMax;
//			if (deltaMax*2 <= realValue*MAX_RATIO){
//				deltaMax *= 2;
//			}
//			else {
//				deltaMax += (realValue*MAX_RATIO-deltaMax+1)/2;
//			}
//		}
//		while(deltaMax - deltaMin > 1){
//			deltaMid = (deltaMax + deltaMin)/2;
//			if (compareErrors(error(realValue + deltaMid, realValue), maxError) <= 0){
//				deltaMin = deltaMid;
//			}
//			else deltaMax = deltaMid;
//		}
//
//		return realValue + deltaMin;
//	}
//
//	@Deprecated
//	public double accumulate (double accumulatedError, int realValue, int observedValue){
//		return accumulate(accumulatedError, error(realValue, observedValue));
//	}
}
