package cgh;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class Simulator {

	private static final double EPSILOIN = 0.0000001;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//		printRatios(10);
		//		System.exit(0);

		Random random = new Random(0);
		int probeNum = 100000;
		double minHitProb = 0.0005;
		double maxHitProb = 0.01;
		double hitProbStep = 0.0005;

		for (int cn = 1; cn<5; ++cn){
			for (int markerNum = 5; markerNum < 100; markerNum += 10){
				System.out.println("Copy number: " + cn + ", markers per segment: " 
						+ markerNum + ", hit probability between " + minHitProb 
						+ " and " + maxHitProb + " with steps of " + hitProbStep);
				for (double hibProb = minHitProb; hibProb <= maxHitProb; hibProb += hitProbStep){
					System.out.print(simulateSegmentStd(markerNum, probeNum, cn, hibProb, random) + ", ");
				}
				
				System.out.println();
				System.out.println();
//					
//					int[] histogram = ratioHistogram(markerNum, probeNum, cn, hibProb, random);
//					double[] split = split(histogram);
//					int minIx = 0;
//					for(; minIx<split.length-1 && split[minIx] <= split[minIx+1]; ++minIx);
//					int maxIx = minIx;
//					for(; maxIx<split.length-1 && split[maxIx+1] == split[maxIx]; ++maxIx);
//					System.out.print("Copy number: " + cn + ", predicted: ");
//					if (maxIx == minIx){
//						System.out.println(minIx);
//					}
//					else{
//						System.out.println(minIx+"-"+maxIx);
//					}
//
//					System.out.println(Arrays.toString(histogram));
//					System.out.println(Arrays.toString(sum(histogram)));
//					//			System.out.println(Arrays.toString(split));
//					System.out.println();
			}
		}
	}

	public static double simulateSegmentStd(int markerNum, int probeNum, int cn, double hibProb, Random random){
		double sum = 0;
		double sumOfSquares = 0;

		for (int i=0; i<markerNum; ++i){
			double ratio = ratioSimulator(probeNum, cn, hibProb, random);
			sum += ratio;
			sumOfSquares += ratio * ratio;
		}

		double average = sum / markerNum;
		return Math.sqrt(sumOfSquares /  markerNum - average*average);
	}

	private static int[] sum(int[] histogram) {
		int[] sum = Arrays.copyOf(histogram, histogram.length);
		for (int i=1; i<sum.length; ++i){
			sum[i] += sum[i-1];
		}
		return sum;
	}

	private static double[] split(int[] histogram) {
		int n = histogram.length;
		double[] split = new double[n+1];
		split[0] = EPSILOIN;
		for (int i = 1; i <= n; ++i){
			split[i] = split[i-1] + histogram[i-1];
		}
		split[n] += EPSILOIN;
		double total = split[n];
		double factor = total * total * 0.25;
		for (int i = 0; i <= n; ++i){
			split[i] = (split[i] * (total-split[i]))/factor;
		}
		return split;
	}

	public static double ratioSimulator(int probes, int cn, double hibProb, Random random){
		double red = 0, green = 0;
		int n = cn+2;
		for (int j = 0; j < probes; ++j){
			if (random.nextDouble() <= hibProb){
				if (random.nextInt(n) < 2){
					++green;
				}
				else{
					++red;
				}
			}
		}
		return red/green*2;
	}

	public static int[] ratioHistogram(int markerNum, int probes, int cn, double hibProb, Random random){
		double[] markers = new double[markerNum];
		for (int i=0; i<markerNum; ++i){
			markers[i] = ratioSimulator(probes, cn, hibProb, random);
		}
		Arrays.sort(markers);
		int last = markerNum;
		for (; last > 0 && Double.isInfinite(markers[last-1]); --last);
		int maxVal = (int) (markers[last-1]+0.5);
		int[] histogram = new int[maxVal+1];
		int ix = 0;
		for (int val = 0; val <= maxVal; ++val){
			for (; ix < last && markers[ix] < val+1; ++ix){
				++histogram[val];
			}
		}
		return histogram;
	}

	public static void printRatios(double n){
		double factor = 4.0/n/n;
		List<Long> ratios = new ArrayList<Long>();
		long currRatio = 0;
		for (int i=1; i<n; ++i){
			long ratio = Math.round(i/(n-i)*2);
			//			System.out.println(i + ", " + (i * (n-i) * factor) + ", " + ratio);
			if (ratio != currRatio){
				ratios.add(ratio);
				currRatio = ratio;
			}
		}
		System.out.print(ratios);
	}


}
