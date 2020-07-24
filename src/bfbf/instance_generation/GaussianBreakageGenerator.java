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

import java.util.Random;

/**
 * A random genome generator who applies normal distribution for choosing 
 * positions of rearrangement events. 
 * 
 * @author Shay Zakov
 *
 */
public class GaussianBreakageGenerator implements BfbBreakageGenerator {

	private Random random;
	private double factor;

	public GaussianBreakageGenerator(Random random, double factor) {
		this.random = random;
		this.factor = factor;
	}

	@Override
	public double randomBreakPosition() {
		double value = 2;
		while (value > 1){
			value = Math.abs(random.nextGaussian()*factor);
		}
		return value;
	}

}
