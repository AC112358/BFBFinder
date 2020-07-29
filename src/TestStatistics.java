import bfbf.*;
import bfbf.weights.ErrorModel;
import bfbf.weights.PoissonErrorModel;
import bfbf.weights.Weights;
import cgh.BFBFileReader;
import cgh.ChromosomeArm;

import java.io.IOException;
import java.util.*;

public class TestStatistics {

    public static void main(String[] args) throws IOException {
        String path = "randomBFBCounts.txt";
        List<int[]> countsLs = BFBFileReader.readSimCounts(path);
        double minWeight = 0.90;
        PoissonErrorModel em = new PoissonErrorModel();
        int pruneMethod = 0;
        int[] candidateCounts = new int[2];
        int[] prunedCandidateCounts = new int[2];


        for (int i = 0; i < countsLs.size(); ++i) {
            int[] counts = countsLs.get(i);
            Weights w = new Weights(counts.length);
//			em.getCountWeights(counts, w.minCounts, w.weights, minWeight);
            w.processWeights();

//			long start = System.currentTimeMillis();

            System.out.print("#" + i + ": ");
            List<Solution1> prunedSolutions = Signature.AllEtaBFBVectors(w, minWeight, 2, prunedCandidateCounts);
            System.out.print("#" + i + ": ");
            List<Solution1> solutions = Signature.AllEtaBFBVectors(w, minWeight, 0, candidateCounts);
            System.out.println(candidateCounts[0] / prunedCandidateCounts[0] + ", " + candidateCounts[1] / prunedCandidateCounts[1]);
            if (solutions.size() != prunedSolutions.size()) {
                System.err.println(solutions.size() + " != " + prunedSolutions.size() + "!!!");
            }
//			System.out.println(solutions.size() + ", time: " + (System.currentTimeMillis() - start));
        }
    }

    public static void testCGP() throws IOException {
        String path = "/BFB/CGP_count_vectors.txt";
//		List<String>[] metaData = new List[2];
//		metaData[0] = new ArrayList<String>();
//		metaData[1] = new ArrayList<String>();
//		List<int[]>[] counts = new List[2];
//		counts[0] = new ArrayList<int[]>();
//		counts[1] = new ArrayList<int[]>();
//		BFBFileReader.readData(metaData, counts,
//				path);

        List<ChromosomeArm> arms = BFBFileReader.readData(path);

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

        Set<int[]> unique = new HashSet<>();

        int[] neighbor;


        int[] minCounts = new int[10];
        double[][] countProbs = new double[10][];
        ErrorModel model = new PoissonErrorModel();

        int minK = 7;
        double minWeight = 0.8;


        for (int i = 0; i < arms.size(); ++i) {
//			for (int j = 0; j<2; ++j){
            int[] currCounts = arms.get(i).getExpectedCounts();
            if (currCounts == null || Solution.effectiveLength(currCounts) < minK) {
                continue;
            }
            //			if (Solution.effectiveLength(arm.getCounts()) < minK){

            ++tested;
            int size = currCounts.length;
            if (minCounts.length < size) {
                minCounts = new int[size];
                countProbs = new double[size][];
            }

            Weights w = new Weights(currCounts.length);
            model.getCountWeights(currCounts, w.minCounts, w.weights, minWeight);
            w.processWeights();

//				model.getCountProbabilities(arm, minCounts, countProbs, minWeightRatio);

//				BFB_Algo.divideByMaxWeight(countWeights);
//				if (arm.isP()){
//					for(int l = size/2 - 1; l >= 0; --l){
//						int tmp = minCounts[l];
//						minCounts[l] = minCounts[size-l-1];
//						minCounts[size-l-1] = tmp;
//
//						double[] temp = countProbs[l];
//						countProbs[l] = countProbs[size-l-1];
//						countProbs[size-l-1] = temp;
//					}
//				}

            boolean hasBFB = false;
            for (ploidy = 1; ploidy < 5; ++ploidy) {
                //			ploidy = arm.estimatePloidy();
                //			arm.setCounts(BFB_Algo.fixPloidy(arm.getCounts(), ploidy));

                int[] fixedMinCounts = BFB_Algo.fixPloidy(minCounts, ploidy);


                //					error = BFB.findMinMaxErrorToValidNeighbor(count, maxError);
                //					neighbor = BFB.minCanberraErrorValidNeighbor(count, maxError*count.length);
                Solution res = BFBCalculator.longestBFBSubstring(fixedMinCounts, countProbs, maxError);
                if (res.getEffectiveLength() > minK) {
                    hasBFB = true;
//						System.out.println(arm);
                    System.out.println(res.toString());
                    int[] effectiveCounts = res.getEffectiveCounts();
                    System.out.println("Effective solution count vector: " + Arrays.toString(effectiveCounts));
                    //					System.out.println(BFBCalculator.searchBFB(effectiveCounts));
                }
                if (hasBFB) {
                    ++matched;
                    System.out.println(arm);
                    System.out.println();
                }
            }
//			}
        }

//		System.out.println("Samples: " + samples);
        System.out.println("Samples with effective length at least " + minK + ": " + tested);
        System.out.println("Samples with solution of effective length at least " + minK + " and error of at most " + maxError + ": " + matched);


    }


}
