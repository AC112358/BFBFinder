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
 * An error model that allows no deviation from the input counts. 
 * 
 * @author zakov
 *
 */
public class NoErrorModel extends ErrorModel {
	
//	@Override
//	public double error(int realValue, int observedValue) {
//		if (observedValue != realValue){
//			return Double.POSITIVE_INFINITY;
//		}
//		else return 0;
//	}
	
	public String toString(){
		return "no error";
	}

//	@Override
//	public double accumulate(double accumulatedError, double newError) {
//		return accumulatedError + newError;
//	}

	@Override
	public int minRealValue(int observedValue, double maxCurrError) {
		return observedValue;
	}

	@Override
	public int maxRealValue(int observedValue, double maxCurrError) {
		return observedValue;
	}

//	@Override
//	public double maxCurrError(double accumulatedError, double maxError) {
//		return 0;
//	}
//
//	@Override
//	public double normlizeError(double error) {
//		// TODO Auto-generated method stub
//		return error;
//	}
//
//	@Override
//	public double deNormlizeError(double normlizedError) {
//		// TODO Auto-generated method stub
//		return normlizedError;
//	}
//
	@Override
	public double weight(int trueCount, int estimatedCount) {
		if (trueCount == estimatedCount){
			return 1;
		}
		else return 0;
	}

//	@Override
//	@Deprecated
//	public void getCountWeights(int[] counts, int[] minCounts,
//			double[][] countProbabilities, double minCountProb) {
//		// TODO Auto-generated method stub
//
//	}

}
