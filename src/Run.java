import bfbf.*;
import bfbf.palindromes.BFBPalindrome;
import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.ErrorModel;
import bfbf.weights.NoErrorModel;
import bfbf.weights.PoissonErrorModel;
import bfbf.weights.Weights;
import cgh.BFBFileReader;
import cgh.CGPFileStatistics;
import cgh.ChromosomeArm;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

//import bfbf.Bfb;

public class Run {

    public static final Random random = new Random();
    public static final int minK = 8;
    public static final int maxK = 3000;//30;
    public static final int permutations = 1;//10000;
    //	public static final double maxError = 0.3;
    private static final NormalDistribution normalDistribution = new NormalDistribution();
    private static final double minWeight = 0.8;
    private static final double weightFactor = 1;
    //	public static final String path = "data/BingellRied.txt"; //"CGP_count_vectors.txt";
    public static String path = "data/Ried/segmentation/";

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        testBignell();
        System.exit(0);


        path = "../BFB/CGP_count_vectors.txt";

        List<ChromosomeArm> arms = BFBFileReader.readData(path);

        //		List<ChromosomeArm> arms = CGHFileHandler.getAllChromosomeArms(path, ".counts.", minWeightRatio, weightFactor);

        //		File dir = new File(path);

        //		path = "data/CGP/CN/NCI-H508.sgm";
        //		ChromosomeArm[][] arms = ChromosomeArm.getArmsFromSgm(path, "NCI-H508", 2, minWeightRatio, weightFactor);

        //		int samples = arms.size();
        int ploidy = 2;

        //		int validCounts = 0;
        int validPermutationAssertedCounts = 0;
        int totalValidPermutations = 0;
        double validPermutationProb;
        int tested = 0, matched = 0;
        int[] count, permutated;
        Collection<int[][]> solutions;
        //		boolean asserted;
        double error;

        SortedSet<double[]> validCounts = new TreeSet<double[]>(new Comparator<double[]>() {

            @Override
            public int compare(double[] arg0, double[] arg1) {
                // comparing p-values:
                if (arg0[3] < arg1[3]) return -1;
                else if (arg0[3] > arg1[3]) return 1;
                    //comparing errors:
                else if (arg0[2] > arg1[2]) return 1;
                else return -1;
            }
        });
        Set<int[]> unique = new HashSet<int[]>();

        int[] neighbor;


        int[] minCounts = new int[10];
        double[][] countProbs = new double[10][];
        ErrorModel model = new PoissonErrorModel();

        List<int[]> allEffectiveCounts = new ArrayList<int[]>();

        for (int i = 0; i < arms.size(); ++i) {
            //			for (int j = 0; j<2; ++j){
            ChromosomeArm arm = arms.get(i); //[i][j];
            int[] expectedCounts = arm.getExpectedCounts();
            if (arm == null || Solution.effectiveLength(expectedCounts) < minK) {
                continue;
            }
            //			if (Solution.effectiveLength(arm.getCounts()) < minK){

            ++tested;
            int size = arm.size();
            if (minCounts.length != size) {
                countProbs = new double[size][];
                Env.returnIntArray(minCounts);
                minCounts = Env.borrowIntArray(size);
            }


            //				BFB_Algo.divideByMaxWeight(countWeights);
            //			if (arm.isP()){
            //				for(int l = size/2 - 1; l >= 0; --l){
            //					int tmp = minCounts[l];
            //					minCounts[l] = minCounts[size-l-1];
            //					minCounts[size-l-1] = tmp;
            //
            //					double[] temp = countProbs[l];
            //					countProbs[l] = countProbs[size-l-1];
            //					countProbs[size-l-1] = temp;
            //				}
            //			}

            boolean hasBFB = false;
            for (ploidy = 1; ploidy < 5; ++ploidy) {
                //			ploidy = arm.estimatePloidy();
                //			arm.setCounts(BFB_Algo.fixPloidy(arm.getCounts(), ploidy));

                int[] fixedExpCounts = BFB_Algo.fixPloidy(expectedCounts, ploidy);
                //				model.getCountProbabilities(fixedExpCounts, minCounts, countProbs, minWeightRatio); ///OLD
                Weights w = model.getWeights(fixedExpCounts, minWeight); //NEW


                //					error = BFB.findMinMaxErrorToValidNeighbor(count, maxError);
                //					neighbor = BFB.minCanberraErrorValidNeighbor(count, maxError*count.length);

                //				Solution res = BFBCalculator.longestBFBSubstring(minCounts, countProbs, 1-minWeightRatio); // OLD
                Solution1 res = Signature.heaviestBFBVector(w, minK, minWeight);


                if (res != null && res.getEffectiveLength() > minK) {
                    if (!hasBFB) {
                        // first ploidy with a solution
                        ++matched;
                        System.out.println();
                        System.out.println(arm);
                    }
                    hasBFB = true;
                    System.out.println("Ploidy: " + ploidy);
                    System.out.println(res.toString());
                    int[] effectiveCounts = res.getEffectiveCounts();
                    allEffectiveCounts.add(effectiveCounts);
                    System.out.println("Effective solution count vector: " + Arrays.toString(effectiveCounts));
                    //					System.out.println(BFBCalculator.searchBFB(effectiveCounts));
                    //					System.out.println(BFB.findAdmittingString(effectiveCounts));
                    //					System.out.println(BFBCalculator.searchBFB(effectiveCounts));

                }
            }
            //			}
        }

        //		System.out.println("Samples: " + samples);
        System.out.println();
        System.out.println("Samples with effective length at least " + minK + ": " + tested);
        System.out.println("Samples with solution of effective length at least " + minK + " and weight at least " + minWeight + ": " + matched);
        System.out.println();
        System.out.println("-------------------------------------------------");
        System.out.println();

        Weights weights = new Weights(0);
        int[] currCandidateCounts = new int[2];
        int[] currPrunedCandidateCounts = new int[2];
        int[] totalCandidateCounts = new int[2];
        int[] totalPrunedCandidateCounts = new int[2];
        int totalSolutions = 0;

        int type = 1;

        long isTime = 0;
        long nonIsTime = 0;
        long start, end;

        model = new NoErrorModel();

        for (int[] counts : allEffectiveCounts) {
            System.out.print(Arrays.toString(counts));

            int size = counts.length;
            if (minCounts.length != size) {
                countProbs = new double[size][];
                Env.returnIntArray(minCounts);
                minCounts = Env.borrowIntArray(size);
            }

            //			model.getCountProbabilities(counts, minCounts, countProbs, minWeight); //OLD
            Weights w = model.getWeights(counts, 1); //NEW
            //			weights.setWeights(minCounts, countProbs); //OLD

            start = System.currentTimeMillis();
            List<BFBPalindrome> palindromes = BFB.allBFBPalindormes(w, 1, w.length());

            //			int size0 = Signature.AllEtaBFBVectors(weights, minWeight, 2, currPrunedCandidateCounts).size(); //OLD
			/*			List<int[]> allHeavyBFBSubVectors = Signature.allHeavyBFBSubVectors(w, minWeight, counts.length); //NEW
			int size0 = allHeavyBFBSubVectors.size(); //NEW
			List<String> allBFBStrings = new ArrayList<>();
			for (int[] currCounts : allHeavyBFBSubVectors){
				w = model.getWeights(currCounts, 1);
				allBFBStrings.addAll(Arrays.asList(Bfb.allBFBStrings(w, 1, w.length())));
			}


			totalSolutions += size0;
			 */
            end = System.currentTimeMillis();
            isTime += end - start;

            //			start = System.currentTimeMillis();
            //			int size1 = Signature.AllEtaBFBVectors(weights, minWeight, 0, currCandidateCounts).size();
            //			if (size0 != size1){
            //				System.err.println(size0 + " != " + size1);
            //			}
            //			end = System.currentTimeMillis();
            //			nonIsTime += end - start;

            totalPrunedCandidateCounts[0] += currPrunedCandidateCounts[0];
            totalPrunedCandidateCounts[1] += currPrunedCandidateCounts[1];
            totalCandidateCounts[0] += currCandidateCounts[0];
            totalCandidateCounts[1] += currCandidateCounts[1];

            //			System.out.print(", vectors: " + size0);
            System.out.println(", strings: " + palindromes.size());
            //			System.out.println(", solutions: " + size0 +
            //					", examined pruned solutions: " + currPrunedCandidateCounts[type] +
            //					", examined un-pruned solutions: " + currCandidateCounts[type] +
            //					", ratio: " + currCandidateCounts[type]/currPrunedCandidateCounts[type]);
        }


        System.out.println("Total solutions: " + totalSolutions +
                ", examined pruned solutions: " + totalPrunedCandidateCounts[type] +
                ", examined un-pruned solutions: " + totalCandidateCounts[type] +
                ", ratio: " + totalCandidateCounts[type] / totalPrunedCandidateCounts[type]);


        System.out.println("Total IS time (sec): " + isTime + ", total non-IS time: "
                + nonIsTime + ", ration: " + nonIsTime / isTime);

        //		System.out.println("Maximum allowed (avarage) error: " + maxError);
        //		System.out.println("Counts with valid neighbors (out of " + tested +" tested samples): " + validCounts.size() + " (unique: " + unique.size() +")");
        //		//		System.out.println("Valid counts for which none of " + permutations + " random permuatations were valid (out of " + tested +" tested samples): " + validPermutationAssertedCounts);
        //		//		System.out.println("Average valid permutations per valid input: " + totalValidPermutations/validCounts);
        //		System.out.println("Name\tChromosome\tOriginal count\tValid neighbor\tError\tPermutation pValue\tNeighbor's string");
        //		for (double[] solution : validCounts){
        //			System.out.print(metaData[0].get((int) solution[0]) + "\t");
        //			System.out.print(metaData[1].get((int) solution[0]) + "\t");
        //			count = counts[(int) solution[1]].get((int) solution[0]);
        //			System.out.print(Arrays.toString(count) + "\t");
        ////			int[] validNeighbor = BFB.minMaxErrorValidNeighbor(count, solution[2]).iterator().next()[0];
        ////			validNeighbor = Arrays.copyOfRange(validNeighbor, 1, validNeighbor.length-1); // trimming added auxiliary entries.
        //			int[] validNeighbor = BFB.minCanberraErrorValidNeighbor(count, maxError*count.length);
        //			System.out.print(Arrays.toString(validNeighbor) + "\t");
        //			System.out.print(solution[2] + "\t");
        //			System.out.print(solution[3] + "\t");
        //					System.out.println(BFB.findAdmittingString(validNeighbor));
        //		}
        //
        //		System.out.println();
        //		System.out.println("Building Markov model...");
        //
        //		int[][] allCounts = merge(counts, minK, maxK);
        //		HMM model = buildHmm(allCounts);
        //
        //		List<int[]> validSimulated = new ArrayList<int[]>();
        //		List<Double> validSimulatedErrors = new ArrayList<Double>();
        //		List<Double> validSimulatedPermutationProbs = new ArrayList<Double>();
        //		unique.clear();
        //
        //		System.out.println("Validating simulated data...");
        //
        //		for (int i=0; i<allCounts.length; ++i){
        //			count = model.sampleWord(allCounts[i].length);
        //			neighbor = BFB.minCanberraErrorValidNeighbor(count, maxError*count.length);
        //			if (neighbor != null){
        //				error = BFB.canberraError(count, neighbor)/count.length;
        //				validSimulated.add(count);
        //				validSimulatedErrors.add(error);
        //				validSimulatedPermutationProbs.add(validPermutationProbability(count, maxError*count.length));
        //				boolean isNew = true;
        //				for (int[] u : unique){
        //					if (Arrays.equals(u, count)){
        //						isNew = false;
        //						break;
        //					}
        //				}
        //				if (isNew) unique.add(count);
        //			}
        //		}
        //
        //		System.out.println("Counts with valid neighbors (out of " + allCounts.length +" tested samples): " + validSimulated.size() + " (unique: " + unique.size() +")");
        //		//		System.out.println("Valid counts for which none of " + permutations + " random permuatations were valid (out of " + tested +" tested samples): " + validPermutationAssertedCounts);
        //		//		System.out.println("Average valid permutations per valid input: " + totalValidPermutations/validCounts);
        //		System.out.println("Original count\tValid neighbor\tError\tPermutation pValue\tNeighbor's string");
        //		for (int i=0; i<validSimulated.size(); ++i){
        //			count = validSimulated.get(i);
        //			System.out.print(Arrays.toString(count) + "\t");
        ////			int[] validNeighbor = BFB.minMaxErrorValidNeighbor(count, validSimulatedErrors.get(i)).iterator().next()[0];
        ////			validNeighbor = Arrays.copyOfRange(validNeighbor, 1, validNeighbor.length-1); // trimming added auxiliary entries.
        //			int[] validNeighbor = BFB.minCanberraErrorValidNeighbor(count, maxError*count.length);//validSimulatedErrors.get(i));
        //			System.out.print(Arrays.toString(validNeighbor) + "\t");
        //			System.out.print(validSimulatedErrors.get(i) + "\t");
        //			System.out.print(validSimulatedPermutationProbs.get(i) + "\t");
        //			System.out.println(BFB.findAdmittingString(validNeighbor));
        //		}
        //
    }


    public static void testBignell2() throws IOException {

        String identifier = "genotypes.csv.gz"; //".genotypes.csv.gz";
        String inputDir = "../BFB/data/CGP/";
        String outputDir = "data/CGP/analyses/CN/";
        File dir = new File(inputDir);
        String[] files = dir.list();
        List<ChromosomeArm> arms = new ArrayList<ChromosomeArm>();

        for (String file : files) {
            if (file.contains(identifier)) {
                String path = inputDir + file;
                System.out.print("Collecting statistics from " + path);
                long start = System.currentTimeMillis();
                CGPFileStatistics statistics = CGPFileStatistics.make(path, "CN");
                statistics.writeToFile(outputDir + file.substring(0, file.indexOf(".") + 1) + "sgm");
                System.out.println(" (" + (System.currentTimeMillis() - start) / 1000 + " sec)");
            }
        }


        String path = "data/CGP/CGP_count_vectors.txt";
        int minK = 10;
        double minWeight = 0.9;
        String outputPath = "data/CGP/analyses/";
        File logFile = new File(outputPath + "picnic_log.txt");
        logFile.createNewFile();
        PrintStream log = new PrintStream(logFile);
        log = System.out;

        List<ChromosomeArm> arms = BFBFileReader.readData(path);
        int tested = 0, matched = 0;

        ErrorModel model = new PoissonErrorModel();
        int maxStrings = 10000;

        ErrorModel noErrors = new NoErrorModel();


        log.println("Input file: " + path);
        log.println("Total number of chromosomal arms: " + arms.size());
        log.println("Error model: " + model.toString());
        log.println("Minimum sub-vector length: " + minK);
        log.println("Minimum sub-vector Weight: " + minWeight);
        log.println("Maximum examined strings per sub-vector: " + maxStrings);

        AllPairwiseBFBStringAccumulator handler = new AllPairwiseBFBStringAccumulator(null, 0, 0, 1, maxStrings);

        PrintStream out = null;

        for (int i = 0; i < arms.size(); ++i) {
            ChromosomeArm arm = arms.get(i);
            int[] expectedCounts = arm.getExpectedCounts();
            if (arm == null || !(arm.getSampleName().contains("H508") && arm.getChromosomeIx() == 14) || Solution.effectiveLength(expectedCounts) < minK) {
                continue;
            }

            ++tested;
            boolean hasBFB = false;

            for (int ploidy = 1; ploidy < 5; ++ploidy) {
                int[] fixedExpCounts = BFB_Algo.fixPloidy(expectedCounts, ploidy);
                Weights w = model.getWeights(fixedExpCounts, minWeight);
                List<int[]> bfbVectors = Signature.allHeavyBFBSubVectors(w, minWeight, minK);
                for (int j = bfbVectors.size() - 1; j >= 0; --j) {
                    if (Solution.effectiveLength(bfbVectors.get(j)) < minK) {
                        bfbVectors.remove(j);
                    }
                }

                if (!bfbVectors.isEmpty()) {
                    handler.clear();

                    if (!hasBFB) {
                        ++matched;
                        log.println();
                        log.println(arm);
                        File outputFile = new File(outputPath + arm.getSampleName()
                                + "_" + arm.getChromosomeIx() + "_" + (arm.isP() ? "P" : "Q") + ".txt");
                        outputFile.createNewFile();
                        out = new PrintStream(outputFile);

                        out = System.err;

                        out.println();
                        out.println(arm);
                        hasBFB = true;
                        handler.ensureLength(w.length());
                    }

                    log.println("Ploidy: " + ploidy + ", number of bfb vectors: " + bfbVectors.size());
                    out.println("Ploidy: " + ploidy + ", number of bfb vectors: " + bfbVectors.size());

                    for (int[] vec : bfbVectors) {
                        out.println(Arrays.toString(vec));
                        int end = Solution.end(vec);
                        handler.set(noErrors.getWeights(vec, 1), Solution.start(vec), end, 1);
                        handler.handle(new PalindromeCollection(), end - 1, 1);
                    }

                    if (handler.getStringCount() >= maxStrings) {
                        log.println("Over " + maxStrings + " bfb strings!");
                        out.println("Over " + maxStrings + " bfb strings!");
                    } else {
                        log.println("Number of bfb strings: " + handler.getStringCount());
                        out.println("Number of bfb strings: " + handler.getStringCount());
                    }
                    out.println();
                    out.println("All pairwise string projections:");
                    handler.setFromTo(0, w.length());
                    handler.printAll(out);

                    log.println();
                    out.println();
                }
            }

            if (hasBFB) {
                out.close();
            }
        }

        log.println();
        log.println("Samples with effective length at least " + minK + ": " + tested);
        log.println("Samples with solution of length at least " + minK + " and weight at least " + minWeight + ": " + matched);
        log.close();
    }


    public static void testBignell() throws IOException {
        String path = "data/CGP/CGP_count_vectors.txt";
        int minK = 8;
        double minWeight = 0.7;
        String outputPath = "data/CGP/analyses2/";
        File logFile = new File(outputPath + "picnic_log.txt");
        logFile.createNewFile();
        PrintStream log = new PrintStream(logFile);
//		log = System.out;

        List<ChromosomeArm> arms = BFBFileReader.readData(path);
        int tested = 0, matched = 0;

        ErrorModel model = new PoissonErrorModel();
        int maxStrings = 10000;

        ErrorModel noErrors = new NoErrorModel();
        int totalBFBs = 0;
        int totalPloidySamples = 0;


        log.println("Input file: " + path);
        log.println("Total number of chromosomal arms: " + arms.size());
        log.println("Error model: " + model.toString());
        log.println("Minimum sub-vector length: " + minK);
        log.println("Minimum sub-vector Weight: " + minWeight);
        log.println("Maximum examined strings per sub-vector: " + maxStrings);

        AllPairwiseBFBStringAccumulator handler =
                new AllPairwiseBFBStringAccumulator(null, 0, 0, 1, maxStrings);

        PrintStream out = null;

        for (int i = 0; i < arms.size(); ++i) {
            ChromosomeArm arm = arms.get(i);
            int[] expectedCounts = arm.getExpectedCounts();
            if (Solution.effectiveLength(expectedCounts) < minK) { // || !(arm.getSampleName().contains("H508") && arm.getChromosomeIx() == 14)
                continue;
            }

            ++tested;
            boolean hasBFB = false;

            for (int ploidy = 1; ploidy < 5; ++ploidy) {
                int[] fixedExpCounts = BFB_Algo.fixPloidy(expectedCounts, ploidy);
                Weights w = model.getWeights(fixedExpCounts, minWeight);
                List<int[]> bfbVectors = Signature.allHeavyBFBSubVectors(w, minWeight, minK);
                for (int j = bfbVectors.size() - 1; j >= 0; --j) {
                    if (Solution.effectiveLength(bfbVectors.get(j)) < minK) {
                        bfbVectors.remove(j);
                    }
                }

                if (!bfbVectors.isEmpty()) {
                    handler.clear();
                    ++totalPloidySamples;
                    totalBFBs += bfbVectors.size();

                    if (!hasBFB) {
                        ++matched;
                        log.println();
                        log.println(arm);
                        File outputFile = new File(outputPath + arm.getSampleName()
                                + "_" + arm.getChromosomeIx() + "_" + (arm.isP() ? "P" : "Q") + ".txt");
                        outputFile.createNewFile();
                        out = new PrintStream(outputFile);

//						out = System.err;

                        out.println();
                        out.println(arm);
                        hasBFB = true;
                        handler.ensureLength(w.length());
                    }

                    log.println("Ploidy: " + ploidy + ", number of bfb vectors: " + bfbVectors.size());
                    out.println("Ploidy: " + ploidy + ", number of bfb vectors: " + bfbVectors.size());

                    for (int[] vec : bfbVectors) {
                        out.println(Arrays.toString(vec));
                        int end = Solution.end(vec);
                        handler.set(noErrors.getWeights(vec, 1), Solution.start(vec), end, 1);
                        handler.handle(new PalindromeCollection(), end - 1, 1);
                    }

                    if (handler.getStringCount() >= maxStrings) {
                        log.println("Over " + maxStrings + " bfb strings!");
                        out.println("Over " + maxStrings + " bfb strings!");
                    } else {
                        log.println("Number of bfb strings: " + handler.getStringCount());
                        out.println("Number of bfb strings: " + handler.getStringCount());
                    }
                    out.println();
                    out.println("All pairwise string projections:");
                    handler.setFromTo(0, w.length());
                    handler.printAll(out);

                    log.println();
                    out.println();
                }
            }

            if (hasBFB) {
                out.close();
            }
        }

        log.println();
        log.println("Samples with effective length at least " + minK + ": " + tested);
        log.println("Samples with solution of length at least " + minK + " and weight at least " + minWeight + ": " + matched);
        log.println("Total ploidy values for such samples: " + totalPloidySamples);
        log.println("Total number of BFB vectors found: " + totalBFBs);
        log.close();
    }


    //	private static HMM buildHmm(int[][] allCounts) {
    //		int maxCount = getMaxCount(allCounts);
    //		int languageSize = (int) (maxCount * (1 + 2 * maxError))+1; //+1 for adding 0 to the language
    //		int n = maxCount+1;
    //		String[] m = new String[languageSize];
    //		double[] emissions = new double[languageSize];
    //
    //		for (int obs = 0; obs<languageSize; ++obs){
    //			m[obs] = ""+obs;
    //		}
    //
    //		HMM model1 = new HMM(m, n);
    //		double std;
    //
    //		for (int state = 0; state < n; ++state){
    //			std = state*maxError;
    //			emissions[0] = normalDistribution.cumulativeProbability(normalDistScale(0.5, state, std));
    //			emissions[languageSize-1] = 1 - normalDistribution.cumulativeProbability(normalDistScale(languageSize-1.5, state, std));
    //
    //			for (int obs=1; obs<languageSize-1; ++obs){
    //				emissions[obs] = normalDistribution.cumulativeProbability(normalDistScale(obs - 0.5, state, std),
    //						normalDistScale(obs + 0.5, state, std));
    //			}
    //
    //			model1.setStateEmissions(state, emissions);
    //
    //		}
    //
    //		model1.baumWelch(allCounts, 1000, 0.01);
    //		return model1;
    //	}

    private static final double normalDistScale(double x, double mu,
                                                double std) {
        return (x - mu) / std;
    }

    private static int getMaxCount(int[][] counts) {
        int maxCount = 0;

        for (int i = 0; i < counts.length; ++i) {
            int[] vec = counts[i];
            for (int j = 0; j < vec.length; ++j) {
                maxCount = Math.max(maxCount, vec[j]);
            }
        }
        return maxCount;
    }

    private static int[][] merge(List<int[]>[] counts, int minLength, int maxLength) {
        int longCountNum = 0;
        for (int i = 0; i < counts[0].size(); ++i) {
            if (counts[0].get(i).length >= minLength && counts[0].get(i).length <= maxLength) {
                ++longCountNum;
            }
            if (counts[1].get(i).length >= minLength && counts[1].get(i).length <= maxLength) {
                ++longCountNum;
            }
        }

        int[][] allCounts = new int[longCountNum][];
        int j = 0;
        for (int i = 0; i < counts[0].size(); ++i) {
            int[] currCount = counts[0].get(i);
            if (currCount.length >= minLength && currCount.length <= maxLength) {
                allCounts[j] = currCount;
                ++j;
            }
            currCount = counts[1].get(i);
            if (currCount.length >= minLength && currCount.length <= maxLength) {
                allCounts[j] = currCount;
                ++j;
            }
        }

        return allCounts;
    }

    //	private static double validPermutationProbability(int[] count, double error) {
    //		int lowerErrorNeighbors;
    //		int[] permutated;
    //		Collection<int[][]> solutions;
    //		lowerErrorNeighbors = 0;
    //		permutated = Arrays.copyOf(count, count.length);
    //		for (int p=0; p<permutations; ++p){
    //			permutate(permutated);
    //			int[] solution = BFB.minCanberraErrorValidNeighbor(permutated, error);
    //			if (solution != null){
    //				++lowerErrorNeighbors;
    //			}
    ////			solutions = BFB.minMaxErrorValidNeighbor(permutated, error);
    ////			if (!solutions.isEmpty()){
    ////				++lowerErrorNeighbors;
    ////			}
    //		}
    //		return ((double)lowerErrorNeighbors)/permutations;
    //	}

    public static void permutate(int[] count) {
        int ix, tmp;

        for (int i = count.length; i > 1; --i) {
            ix = random.nextInt(i);
            tmp = count[i - 1];
            count[i - 1] = count[ix];
            count[ix] = tmp;
        }
    }

}
