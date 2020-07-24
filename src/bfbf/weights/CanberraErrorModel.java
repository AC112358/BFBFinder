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
 * The distance between an observed value q and a real value r is defined by
 * |r - q| / (r + q). The distance between two vectors is defined to be the
 * sum of entry-wise distances. 
 * 
 * @author Shay Zakov
 *
 */
public class CanberraErrorModel extends ErrorModel {

//	@Override
//	public double error(int realValue, int observedValue) {
//		return Math.abs(realValue - observedValue) /
//				((double) realValue + observedValue);
//	}
//
//	@Override
//	public double accumulate(double accumulatedError, double newError) {
//		return accumulatedError + newError;
//	}

	@Override
	public double weight(int trueCount, int estimatedCount) {
		return 1 - Math.abs(trueCount - estimatedCount) /
				((double) trueCount + estimatedCount);
	}

	@Override
	public int minRealValue(int observedValue, double maxCurrError) {
		return (int) Math.max(1, 
				Math.ceil(observedValue*(1-maxCurrError)/(1+maxCurrError)));
	}

	@Override
	public int maxRealValue(int observedValue, double maxCurrError) {
		return (int) Math.min(3*observedValue, 
				Math.floor(observedValue*(1+maxCurrError)/(1-maxCurrError)));
	}

//	@Override
//	public double maxCurrError(double accumulatedError, double maxError) {
//		return maxError - accumulatedError;
//	}
//
//	@Override
//	public double normlizeError(double error) {
//		return error;
//	}
//
//	@Override
//	public double deNormlizeError(double normlizedError) {
//		return normlizedError;
//	}

}
