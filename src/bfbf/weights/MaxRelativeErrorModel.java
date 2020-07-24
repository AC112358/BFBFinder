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
 * For a maximum relative error {@code m}, the ratio between a pertubated value
 * {@code q} and the real value {@code r} satisfies {@code 1/(1+m) <= q/r <= 1+m}.
 * 
 * 
 * 
 * @author Shay Zakov
 *
 */
public class MaxRelativeErrorModel extends ErrorModel {
	
//	@Override
//	public double error(int realValue, int observedValue) {
//		if (observedValue <= realValue){
//			return ((double) realValue)/observedValue - 1;
//		}
//		else return ((double) observedValue)/realValue - 1;
//	}

//	@Override
//	public double accumulate(double accumulatedError, double newError) {
//		return Math.max(accumulatedError, newError);
//	}

	@Override
	public double weight(int trueCount, int estimatedCount) {
		if (estimatedCount >= trueCount){
			return ((double) trueCount)/estimatedCount;
		}
		else return ((double) estimatedCount)/trueCount;
	}

	@Override
	public int minRealValue(int observedValue, double maxCurrError) {
		return (int) Math.ceil((observedValue/(1+maxCurrError)));
	}

	@Override
	public int maxRealValue(int observedValue, double maxCurrError) {
		return (int) Math.floor((observedValue*(1+maxCurrError)));
	}

//	@Override
//	public double maxCurrError(double accumulatedError, double maxError) {
//		return maxError;
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
