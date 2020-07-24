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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * A randome genome generator.
 * 
 * @author Shay Zakov
 *
 */
public class RandomSampleGenerator {

	private static final int BFB_CYCLE = 0;
	private static final int DUPLICATION = 1;
	private static final int TRANSLOCATION = 2;
	private static final int DELETION = 3;
	private static final int INVERSION = 4;

	private static final int BFB_START = 0;
	private static final int BFB_END = 1;
	private static final int BFB_CONSECUTIVE = 2; 
	private static final int BFB_SCATTERED = 3;

	private static PrintStream outStream = System.out; 
	private static PrintStream logSteam = null; 

	private static long genomeLength = 100000000;
	private static long minSegmentLength = 5000;
	private static int iters = 5;
	private static String bfbLocationStr = "start";
	private static int bfbLocation;
	private static int bfbCycles = 6;
	private static String bfbBreakageDistribution = "gaussian"; 
	private static double gaussianFactor = 0.1;
	private static int otherRearrangements = 5;
	private static double duplicationWeight = 1;
	private static double translocationWeight = 1;
	private static double deletionWeight = 1;
	private static double inversionWeight = 1;
	private static double dupInvertionProb = 0.5;
	private static double tandemDupProb = 0.5;
	private static int ploidy = 2;

	private static long segmentLengthMean = 5000; 
	private static long segmentVar = 1000; 

	private static Long randomSeed = null;
	private static double duplicationConst;
	private static double translocationConst;
	private static double deletionConst;
	private static Random random; 
	private static BfbBreakageGenerator bfbBreakageGenerator;

	private static final int MAX_ATTEMPTS = 100; 

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		setEnvironment(args); 
		printEnvironment(outStream);
		if (logSteam != null){
			printEnvironment(logSteam);
		}

		SegmentedChromosome genome = new SegmentedChromosome(genomeLength, ploidy);
		List<Integer> operations = new ArrayList<Integer>();


		for (int i=0; i<iters; ++i){
			log(logSteam, "\n---------- Iteration " + i + " ----------------------");
			System.out.println(i);

			generateBreakpoints(genome, operations, bfbLocation, bfbCycles, otherRearrangements, 
					duplicationConst, deletionConst, random, bfbBreakageGenerator);


			outStream.println(Arrays.toString(genome.getPloidSegments()));
			outStream.println(Arrays.toString(genome.getCounts()));
			outStream.println();
		}
	}


	public static void generateBreakpoints(SegmentedChromosome genome,
			List<Integer> operations, int bfbLocation, int bfbCycles, int otherRearrangements,
			double duplicationConst, double deletionConst, Random random, BfbBreakageGenerator bfbBreakageGenerator) {
		boolean isInverted;
		boolean isTandem;

		genome.clear();
		setOperationSequence(operations, bfbLocation, bfbCycles, otherRearrangements, 
				duplicationConst, deletionConst, random);

		for (int j=0; j<operations.size(); ++j){
			String operationType = null;
			long segmentLength;
			int currOperation = operations.get(j);
			if (currOperation == BFB_CYCLE){
				segmentLength = genome.getBreakpoint(0, bfbBreakageGenerator.randomBreakPosition());
				genome.bfbCycle(0, segmentLength);
				operationType = "BFB cycle on chromosome 0, duplicated suffix length " + segmentLength;
			}

			else{
				int ploidIx = random.nextInt(ploidy);
				int attempt = 0;
				do{
					segmentLength = (long) (random.nextGaussian() * segmentVar + segmentLengthMean);
					++attempt;
				} while ((segmentLength >= genome.ploidLength(ploidIx) || segmentLength <= 0) 
						&& attempt < MAX_ATTEMPTS);
				if (segmentLength >= genome.ploidLength(ploidIx) || segmentLength <= 0){
					continue;
				}
				double relBreakpoint = random.nextDouble();
				long breakpoint = genome.getBreakpoint(ploidIx, relBreakpoint, segmentLength);
				long insertionPoint = genome.getBreakpoint(ploidIx, relBreakpoint, segmentLength);
				if (insertionPoint >= breakpoint){
					insertionPoint += segmentLength;
				}


				switch (currOperation){
				case DUPLICATION:
					isInverted = random.nextDouble() < dupInvertionProb;
					isTandem = random.nextDouble() < tandemDupProb;

					if (isInverted){
						if (isTandem){
							operationType = "Inverted tandem duplication on chromosome " 
									+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength;
						}
						else{
							operationType = "Inverted duplication on chromosome " 
									+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength
									+ ", insertion position " + insertionPoint;
						}
					}
					else if (isTandem){
						operationType = "Tandem duplication on chromosome " 
								+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength;
					}
					else {
						operationType = "Duplication on chromosome " 
								+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength
								+ ", insertion position " + insertionPoint;
					}

					if (isTandem){
						genome.tandemDuplicate(ploidIx, breakpoint, segmentLength, isInverted);
					}
					else{
						genome.duplicate(ploidIx, breakpoint, segmentLength, insertionPoint, isInverted);
					}
					break;

				case TRANSLOCATION:
					isInverted = random.nextDouble() < dupInvertionProb;

					if (isInverted){
						operationType = "Inverted translocation on chromosome " 
								+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength
								+ ", insertion position " + insertionPoint;
					}
					else {
						operationType = "Translocation on chromosome " 
								+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength
								+ ", insertion position " + insertionPoint;
					}

					genome.translocate(ploidIx, breakpoint, segmentLength, insertionPoint, isInverted);
					break;

				case DELETION:
					operationType = "Deletion on chromosome " 
							+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength;
					genome.delete(ploidIx, breakpoint, segmentLength);
					break;

				case INVERSION:
					operationType = "Inversion on chromosome " 
							+ ploidIx +", starting at " + breakpoint + ", segment length " + segmentLength;
					genome.invert(ploidIx, breakpoint, segmentLength);
					break;

				default:
					throw new RuntimeException("Illegal operation.");
				}
			}
			log(logSteam, operationType + "\n"+ genome);

		}

		if (genome.unifyConsecutiveDeletedSegments()){
			log(logSteam, "Unifying consecutive deleted segments\n" + genome);
		}
	}


	private static void setEnvironment(String[] args) {
		Map<String, String> configuration = readConfiguration(args);

		String value = configuration.get("CONFIG");
		if (value != null){
			readConfigFile(value, configuration);
			readConfiguration(args, configuration);
		}

		value = configuration.get("OUTPUT");
		if (value != null){
			try {
				outStream = new PrintStream(value);
			} catch (FileNotFoundException e) {
				System.err.println("Cannot write to " + value);
			}
		}

		value = configuration.get("LOG");
		if (value != null){
			try {
				logSteam = new PrintStream(value);
			} catch (FileNotFoundException e) {
				System.err.println("Cannot write to " + value);
			}
		}


		random = new Random();
		if (randomSeed == null){
			randomSeed = random.nextLong();
		}
		random.setSeed(randomSeed);

		value = configuration.get("GENOME_LENGTH");
		if (value != null){
			genomeLength = Long.parseLong(value);
		}

		value = configuration.get("PLOIDY");
		if (value != null){
			ploidy = Integer.parseInt(value);
		}

		value = configuration.get("MIN_SEGMENT_LENGTH");
		if (value != null){
			minSegmentLength = Long.parseLong(value);
		}

		value = configuration.get("ITERATIONS");
		if (value != null){
			iters = Integer.parseInt(value);
		}

		value = configuration.get("BFB_LOCATION");
		if (value != null){
			if (value.equals("start")){
				bfbLocation = BFB_START;
				bfbLocationStr = value;
			}
			else if (value.equals("end")){
				bfbLocation = BFB_END;
				bfbLocationStr = value;
			}
			else if (value.equals("consecutive")){
				bfbLocation = BFB_CONSECUTIVE;
				bfbLocationStr = value;
			}
			else if (value.equals("scattered")){
				bfbLocation = BFB_SCATTERED;
				bfbLocationStr = value;
			}
		}

		value = configuration.get("BFB_CYCLES");
		if (value != null){
			bfbCycles = Integer.parseInt(value);
		}

		value = configuration.get("GAUSSIAN_FACTOR");
		if (value != null){
			gaussianFactor = Double.parseDouble(value);
		}

		value = configuration.get("BFB_BREAKAGE_DISTRIBUTION");
		if (value != null){
			if (value.equals("gaussian")){
				bfbBreakageDistribution = value;
				bfbBreakageGenerator = new GaussianBreakageGenerator(random, gaussianFactor);
			}
			else{
				bfbBreakageDistribution = "uniform";
				bfbBreakageGenerator = new UniformBreakageGenerator(random);
			}
		}


		value = configuration.get("OTHER_REARRANGEMENTS");
		if (value != null){
			otherRearrangements = Integer.parseInt(value);
		}

		value = configuration.get("DUPLICATION_WEIGHT");
		if (value != null){
			duplicationWeight = Double.parseDouble(value);
		}

		value = configuration.get("TRANSLOCATION_WEIGHT");
		if (value != null){
			translocationWeight = Double.parseDouble(value);
		}

		value = configuration.get("DELETION_WEIGHT");
		if (value != null){
			deletionWeight = Double.parseDouble(value);
		}

		value = configuration.get("INVERSION_WEIGHT");
		if (value != null){
			inversionWeight = Double.parseDouble(value);
		}

		value = configuration.get("TANDEM_DUP_PROB");
		if (value != null){
			tandemDupProb = Double.parseDouble(value);
		}

		value = configuration.get("DUP_INVERSION_PROB");
		if (value != null){
			dupInvertionProb = Double.parseDouble(value);
		}

		value = configuration.get("REARRANGED_SEGMENT_MEAN_LENGTH");
		if (value != null){
			segmentLengthMean = Long.parseLong(value);
		}

		value = configuration.get("REARRANGED_SEGMENT_VAR");
		if (value != null){
			segmentVar = Long.parseLong(value);
		}

		value = configuration.get("RANDOM_SEED");
		if (value != null){
			randomSeed = Long.parseLong(value);
		}

		double totalWeight = duplicationWeight + translocationWeight + deletionWeight + inversionWeight;
		duplicationConst = duplicationWeight/totalWeight;
		translocationConst = duplicationConst + translocationWeight/totalWeight;
		deletionConst = translocationConst + deletionWeight/totalWeight;
	}

	private static void printEnvironment(PrintStream outStream) {
		outStream.println("# ITERATIONS:" + iters);
		outStream.println("# GENOME_LENGTH:" + genomeLength);
		outStream.println("# PLOIDY:" + ploidy);
		outStream.println("# MIN_SEGMENT_LENGTH:" + minSegmentLength);
		outStream.println("# BFB_LOCATION:" +  bfbLocationStr);
		outStream.println("# BFB_CYCLES:" + bfbCycles);
		outStream.println("# BFB_BREAKAGE_DISTRIBUTION:" + bfbBreakageDistribution);
		if (bfbBreakageDistribution.equals("gaussian")){
			outStream.println("# GAUSSIAN_FACTOR:" + gaussianFactor);
		}
		outStream.println("# OTHER_REARRANGEMENTS:" + otherRearrangements);
		outStream.println("# DUPLICATION_WEIGHT:" + duplicationWeight);
		outStream.println("# TRANSLOCATION_WEIGHT:" + translocationWeight);
		outStream.println("# DELETION_WEIGHT:" + deletionWeight);
		outStream.println("# INVERSION_WEIGHT:" + inversionWeight);
		outStream.println("# duplication probability:" + duplicationConst);
		outStream.println("# translocation probability:" + (translocationConst - duplicationConst));
		outStream.println("# deletion probability: " + (deletionConst - translocationConst));
		outStream.println("# inversion probability: " + (1 - deletionConst));
		outStream.println("# TANDEM_DUP_PROB:" + tandemDupProb);
		outStream.println("# DUP_INVERSION_PROB:" + dupInvertionProb);
		outStream.println("# REARRANGED_SEGMENT_MEAN_LENGTH:" + segmentLengthMean);
		outStream.println("# REARRANGED_SEGMENT_VAR:" + segmentVar);

		outStream.println("# RANDOM_SEED:" + randomSeed);
		outStream.println("# ------------------------------------------------------------\n");
	}


	private static void log(PrintStream logSteam, String string) {
		if (logSteam != null){
			logSteam.println(string);
		}
	}

	private static void readConfiguration(String[] args, Map<String, String> config) {
		for (String str: args){
			String[] pair = str.trim().split(":");
			config.put(pair[0], pair[1]);
		}
	}

	private static Map<String, String> readConfiguration(String[] args) {
		Map<String, String> config = new HashMap<String, String>();
		readConfiguration(args, config);
		return config;
	}

	private static void readConfigFile(String path, Map<String, String> config) {
		File file = new File(path);
		if (file.isFile()){
			try {
				@SuppressWarnings("resource")
				BufferedReader br = new BufferedReader(new FileReader(file));
				for (String line = br.readLine(); line !=null; line = br.readLine()){
					String[] pair = line.trim().split(":");
					config.put(pair[0], pair[1]);
				}
			} 
			catch (Exception e) {
				System.err.println("Error when trying to read " 
						+ file.toString() + ":\n" + e.getMessage()); 
			}
		}
	}

	private static void setOperationSequence(List<Integer> operations,
			int bfbLocation, int bfbCycles, int otherRearrangements,
			double duplicationConst, double deletionConst, Random random) {

		operations.clear();
		int currOperation;
		double randomValue;
		for (int i=0; i<otherRearrangements; ++i){
			randomValue = random.nextDouble();
			if (randomValue < duplicationConst){
				currOperation = DUPLICATION;
			}
			else if (randomValue < translocationConst){
				currOperation = TRANSLOCATION;
			}
			else if (randomValue < deletionConst){
				currOperation = DELETION;
			}
			else currOperation = INVERSION;

			operations.add(currOperation);
		}

		int bfbIx;
		if (bfbLocation == BFB_START){
			bfbIx = 0;
		}
		else if (bfbLocation == BFB_CONSECUTIVE){
			bfbIx = (int) (otherRearrangements*random.nextDouble());
		}
		else bfbIx = otherRearrangements;

		for (int i=0; i<bfbCycles; ++i){
			operations.add(bfbIx, BFB_CYCLE);
		}

		if (bfbLocation == BFB_SCATTERED){
			Collections.shuffle(operations);
		}
	}

//	public static int[] addNoize(int[] counts, double noize){
//		for (int i=1; i<counts.length; ++i){
//			double cdf = Math.random()*noize*2 - noize + 0.5;
//			double sum = 0;
//			int j=1;
//			for (; sum < cdf; ++j){
//				sum += Math.exp(PoissonErrorModel.poissonProbabilityApproximation(counts[i], j));
//			}
//			counts[i] = j-1;
//		}
//		return counts;
//	}

}
