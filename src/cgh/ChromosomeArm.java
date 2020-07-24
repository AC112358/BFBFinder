package cgh;

import bfbf.Env;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChromosomeArm {
	
	
	protected String sampleName;
	protected int chromosomeIx;
	protected boolean isP;
//	protected int[] counts;
//	protected int[] minCounts;
//	protected double[][] countWeights;
	protected int[] starts, ends, probes;
	protected double[] logSignalMeans, logSignalVars; 
	
	
	public ChromosomeArm(String sampleName, int chromosomeIx, boolean isP,
			int[] counts) {
		this.sampleName = sampleName;
		this.chromosomeIx = chromosomeIx;
		this.isP = isP;
		this.starts = new int[counts.length];
		this.ends = new int[counts.length];
		this.probes = new int[counts.length];
		this.logSignalMeans = new double[counts.length];
		this.logSignalVars = new double[counts.length];
		Arrays.fill(probes, 1);
		
		for (int i=0; i<counts.length; ++i){
			starts[i] = i;
			ends[i] = i;
			logSignalMeans[i] = Math.log(counts[i]);
			logSignalVars[i] = logSignalMeans[i];
		}
	}
	
	public ChromosomeArm(String sampleName, int chromosomeIx, boolean isP,
			int[] starts, int[] ends, int[] probes, double[] logSignalMeans,
			double[] logSignalVars) {
		super();
		this.sampleName = sampleName;
		this.chromosomeIx = chromosomeIx;
		this.isP = isP;
		this.starts = starts;
		this.ends = ends;
		this.probes = probes;
		this.logSignalMeans = logSignalMeans;
		this.logSignalVars = logSignalVars;
	}
	
	public ChromosomeArm(String sampleName, int chromosomeIx, boolean isP,
			List<Integer> startLs, List<Integer> endLs, List<Integer> probeLs, 
			List<Double> logRatioMeanLs, List<Double> logRatioSdLs) {
		this.sampleName = sampleName;
		this.chromosomeIx = chromosomeIx;
		this.isP = isP;

		int size = startLs.size();
		starts = new int[size]; 
		ends = new int[size];
		probes = new int[size];
		logSignalMeans = new double[size];
		logSignalVars = new double[size];
		
		for (int i=0; i<size; ++i){
			starts[i] = startLs.get(i); 
			ends[i] = endLs.get(i);
			probes[i] = probeLs.get(i);
			logSignalMeans[i] = logRatioMeanLs.get(i);
			logSignalVars[i] = logRatioSdLs.get(i);
		}

	}


//	public ChromosomeArm(String sampleName, int chromosomeIx, boolean isP, int[] counts) {
//		super();
//		this.sampleName = sampleName;
//		this.chromosomeIx = chromosomeIx;
//		this.isP = isP;
//		this.counts = counts;
//	}
	
	public String getSampleName() {
		return sampleName;
	}
	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}
	public int getChromosomeIx() {
		return chromosomeIx;
	}
	public void setChromosomeIx(int chromosomeIx) {
		this.chromosomeIx = chromosomeIx;
	}
	
	public int[] getExpectedCounts() {
		int[] counts = new int[logSignalMeans.length];
		for (int i=0; i<logSignalMeans.length; ++i){
			int count = (int) (Math.round(Math.exp(logSignalMeans[i])));
			if (logSignalMeans[i] > Math.log(count * (count+1))){
				counts[i] = count+1;
			}
			else{
				counts[i] = count;
			}
		}
		return counts;
	}
	
	public void getCountWeights(TIntList minCounts, List<TDoubleList> coutWeights, double minWeight, double beta){
		int[] counts = getExpectedCounts();
		minCounts.clear();
		TDoubleList currWeights;
		int count;
		double weight, sd, diff, weightFactor;
		
		for (int i=0; i<logSignalMeans.length; ++i){
			if (coutWeights.size() == i){
				currWeights = new TDoubleArrayList();
				coutWeights.add(currWeights);
			}
			else {
				currWeights = coutWeights.get(i);
				currWeights.clear();
			}

			count = counts[i]; 
			diff = Math.log(count) - logSignalMeans[i];
			sd = Math.sqrt(logSignalVars[i] + diff*diff);
			weightFactor = 1 / Math.exp(-sd *  beta);
			weight = 1;
			while (weight >= minWeight){
				currWeights.add(weight);
				--count;
				diff = Math.log(count) - logSignalMeans[i];
				sd = Math.sqrt(logSignalVars[i] + diff*diff);
				weight = Math.exp(-sd *  beta) * weightFactor;
			}

			minCounts.add(counts[i] - currWeights.size() + 1);
			currWeights.reverse();

			count = counts[i]+1;
			diff = Math.log(count) - logSignalMeans[i];
			sd = Math.sqrt(logSignalVars[i] + diff*diff);
			weight = Math.exp(-sd *  beta) * weightFactor;
			while (weight >= minWeight){
				currWeights.add(weight);
				++count;
				diff = Math.log(count) - logSignalMeans[i];
				sd = Math.sqrt(logSignalVars[i] + diff*diff);
				weight = Math.exp(-sd *  beta) * weightFactor;
			}
		}
	}
	
	
	public boolean isP() {
		return isP;
	}
	
	public void setIsP(boolean isP) {
		this.isP = isP;
	}
	
	public String toString(){
		return sampleName + " " + chromosomeIx + (isP ? " p\n" : " q\n") + Arrays.toString(getExpectedCounts());
	}

//	public double[][] getCountProbabilities() {
//		return countWeights;
//	}

//	public void setCountProbabilities(double[][] countProbabilities) {
//		this.countWeights = countProbabilities;
//	}

	public int[] getStart() {
		return starts;
	}

	public void setStart(int[] start) {
		this.starts = start;
	}

	public int[] getEnd() {
		return ends;
	}

	public void setEnd(int[] end) {
		this.ends = end;
	}
	
	

//	public void setCounts(List<Integer> minCountLs,
//			List<List<Double>> weights, List<Integer> startLs,
//			List<Integer> endLs) {
//		int size = minCountLs.size();
//		counts = new int[size];
//		minCounts = new int[size];
//		countWeights = new double[size][]; 
//		starts = new int[size]; 
//		ends = new int[size];
//		for (int i=0; i<size; ++i){
//			minCounts[i] = minCountLs.get(i);
//			List<Double> currWeights = weights.get(i);
//			int currSize = currWeights.size();
//			countWeights[i] = new double[currSize];
//			
//			int bestProbIx = 0;
//			
////			int[] ranks = Env.borrowIntArray(currSize);
////			Env.rank(currProbabilities, ranks);
//			
//			double sum = 0;
//
//			for (int j=0; j<currSize; ++j){
//				countWeights[i][j] = currWeights.get(j);
//				sum += countWeights[i][j];
//				if (countWeights[i][bestProbIx] < countWeights[i][j]){
//					bestProbIx = j;
//				}
//			}
//
//			for (int j=0; j<currSize; ++j){
//				countWeights[i][j] /= sum;
//			}
//			
//			counts[i] = minCounts[i] + bestProbIx; 
//			
////			Env.returnIntArray(ranks);
//			starts[i] = startLs.get(i); 
//			ends[i] = endLs.get(i);
//		}
//	}
	
//	public int estimatePloidy(){
//		int maxPossibleCount = getMaxPossibleCount();
//		double[] ploidyWeights = new double[maxPossibleCount+1];
//		for (int i=0; i<minCounts.length; ++i){
//			int length = ends[i] - starts[i];
//			for (int j=0; j<countWeights[i].length; ++j){
//				ploidyWeights[minCounts[i]+j] += length*countWeights[i][j];
//			}
//		}
//		
//		int ploidy = 0;
//		for (int i=1; i<= maxPossibleCount; ++i){
//			if (ploidyWeights[ploidy] < ploidyWeights[i]){
//				ploidy = i;
//			}
//		}		
//		return ploidy;
//	}

//	public int getMaxPossibleCount() {
//		int maxPossibleCount = 0;
//		for (int i=0; i<minCounts.length; ++i){
//			maxPossibleCount = Math.max(maxPossibleCount, minCounts[i] + countWeights[i].length-1);
//		}
//		return maxPossibleCount;
//	}

//	public int[] getMinCounts() {
//		return minCounts;
//	}
//
//	public void setMinCounts(int[] minCounts) {
//		this.minCounts = minCounts;
//	}
//
//	public double[][] getCountWeights() {
//		return countWeights;
//	}
//
//	public void setCountWeights(double[][] countWeights) {
//		this.countWeights = countWeights;
//	}
//
//	public void setP(boolean isP) {
//		this.isP = isP;
//	}
	
	public static ChromosomeArm[][] getArmsFromSgm(String sgmFile, String sampleName, 
			double logBase, double minWeightRation, double weightFactor) throws IOException{
		ChromosomeArm[][] arms = new ChromosomeArm[25][2];

		BufferedReader br = Env.getBufferedReader(sgmFile);

		List<Integer> startLs = new ArrayList<Integer>();
		List<Integer> endLs = new ArrayList<Integer>();
		List<Integer> probeLs = new ArrayList<Integer>();
		List<Double> logRatioMeanLs = new ArrayList<Double>();
		List<Double> logRatioSdLs = new ArrayList<Double>();
		
		double logFactor = Math.log(logBase);

		
//		String gIntP = "([\\d]+)\\t";
//		String intP = "[\\d]+\\t";
//		String gDoubleP = "("+doublePatternStr+")\\t";
//
//		Pattern p = Pattern.compile(intP + gIntP + gIntP + gIntP + intP + gIntP + gDoubleP + gDoubleP + gDoubleP + gDoubleP);

		
		Pattern p = Pattern.compile("\\S+\\s+\\S+\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+");

		int chrGroup = 1;
		int startGroup = 2;
		int endGroup = 3;
		int signalsGroup = 4;
		int meanGroup = 5;
		int sdGroup = 6;

		Matcher matcher;

		String chrStr = null;
		int chrIx, segStart, segEnd, signals, countMedian;
		double mean, maxWeight, currLog, logVar, sd, std;

		int currChromosome = -1;
		boolean isP = true;
		ChromosomeArm currArm = null;

		String line = br.readLine();
		line = br.readLine();
		while (line != null){
			line += "\t";
			matcher = p.matcher(line);
			if (matcher.find()){
				chrStr = matcher.group(chrGroup);
				chrIx = Integer.parseInt(chrStr);
				segStart = Integer.parseInt(matcher.group(startGroup));
				segEnd = Integer.parseInt(matcher.group(endGroup));
				signals = Integer.parseInt(matcher.group(signalsGroup));
				mean = logFactor * Double.parseDouble(matcher.group(meanGroup));
				sd = logFactor * Double.parseDouble(matcher.group(sdGroup)); 
//				sd *= sd;
				
//				countMedian = (int) Math.exp(mean);
//				double lowerLog = Math.log(countMedian);
//				double upperLog = Math.log(countMedian+1);
//				if (mean - lowerLog < upperLog - mean){
//					currLog = lowerLog;
//				}
//				else{
//					currLog = upperLog;
//					++countMedian;
//				}
//
//				logVar = sd + (mean-currLog) * (mean-currLog);
//				std = Math.sqrt(logVar);
//				maxWeight = Math.exp(-std*weightFactor);
//				if (maxWeight == 0){
//					// segment is ignored due to high variance
//					line = br.readLine();
//					continue;
//				}

//				double currWeight = maxWeight;
//				List<Double> currProbs = new ArrayList<Double>();
//
//				for (int i = countMedian-1; i >= -1 && currWeight/maxWeight >= minWeightRation; --i){
//					currProbs.add(currWeight);
//					currLog = Math.log(i);
//					logVar = sd + (mean-currLog) * (mean-currLog);
//					std = Math.sqrt(logVar);
//					currWeight = Math.exp(-std*weightFactor);
//				}
//				
//				int minCount = countMedian - currProbs.size() + 1;
//				Collections.reverse(currProbs);
//				currLog = Math.log(countMedian+1);
//				logVar = sd + (mean-currLog) * (mean-currLog);
//				std = Math.sqrt(logVar);
//				currWeight = Math.exp(-std*weightFactor);
//				for (int i = countMedian+2; currWeight/maxWeight >= minWeightRation; ++i){
//					currProbs.add(currWeight);
//					currLog = Math.log(i);
//					logVar = sd + (mean-currLog) * (mean-currLog);
//					std = Math.sqrt(logVar);
//					currWeight = Math.exp(-std*weightFactor);
//				}

				if (currChromosome != chrIx || (isP && CGHFileHandler.centromerPositions[chrIx] < segEnd)){
					if (currChromosome == chrIx){
						if(CGHFileHandler.centromerPositions[chrIx] > segStart){
							// current segment spans the centromer
							startLs.add(segStart);
							endLs.add(CGHFileHandler.centromerPositions[chrIx]);
							probeLs.add((int) (signals * ((double)CGHFileHandler.centromerPositions[chrIx] -  segStart) / (segEnd - segStart)));
							logRatioMeanLs.add(mean);
							logRatioSdLs.add(sd);
						}

						// finalizing the p-arm of this chromosome:
						currArm = new ChromosomeArm(sampleName, currChromosome, true, startLs, endLs, probeLs, logRatioMeanLs, logRatioSdLs);
						arms[chrIx][0] = currArm;
						isP = false;
					}
					else{
						// finalizing the arm of the previous chromosome:
						if (currChromosome >= 0){
							currArm = new ChromosomeArm(sampleName, currChromosome, isP, startLs, endLs, probeLs, logRatioMeanLs, logRatioSdLs);
							if (isP){
								arms[currChromosome][0] = currArm;
							}
							else{
								arms[currChromosome][1] = currArm;
							}
						}
						currChromosome = chrIx;
						isP = CGHFileHandler.centromerPositions[chrIx] > segStart;
						
						// checking if the new segment spans the centromer of the new chromosome:
						if (isP && CGHFileHandler.centromerPositions[chrIx] < segEnd){
							startLs.add(segStart);
							endLs.add(CGHFileHandler.centromerPositions[chrIx]);
							probeLs.add((int) (signals * ((double)CGHFileHandler.centromerPositions[chrIx] -  segStart) / (segEnd - segStart)));
							logRatioMeanLs.add(mean);
							logRatioSdLs.add(sd);
							currArm = new ChromosomeArm(sampleName, currChromosome, true, startLs, endLs, probeLs, logRatioMeanLs, logRatioSdLs);
							arms[chrIx][0] = currArm;
							isP = false;
						}
					}

					// clearing the lists for the next arm  
					startLs.clear();
					endLs.clear();
					probeLs.clear();
					logRatioMeanLs.clear();
					logRatioSdLs.clear();
				}
				
				if (!isP && CGHFileHandler.centromerPositions[chrIx] > segStart){
					startLs.add(CGHFileHandler.centromerPositions[chrIx]);
					probeLs.add((int) (signals * (segEnd - (double)CGHFileHandler.centromerPositions[chrIx]) / (segEnd - segStart)));
				}
				else{
					startLs.add(segStart);
					probeLs.add(signals);
				}
				endLs.add(segEnd);
				logRatioMeanLs.add(mean);
				logRatioSdLs.add(sd);
			}
			line = br.readLine();
		}

		currArm = new ChromosomeArm(sampleName, currChromosome, isP, startLs, endLs, probeLs, logRatioMeanLs, logRatioSdLs);
		if (isP){
			arms[currChromosome][0] = currArm;
		}
		else{
			arms[currChromosome][1] = currArm;
		}

		return arms;
	}

	public int size() {
		return starts.length;
	}


}
