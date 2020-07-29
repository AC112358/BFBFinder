package cgh;

import bfbf.Env;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class CGHFileHandler {

    public static final char CSV_SEPARATOR = '\t';

    public static final int[] centromerPositions = {
            0,
            124300000, // chr1
            93300000, // chr2
            91700000, // chr3
            50700000, // chr4
            47700000, // chr5
            60500000, // chr6
            59100000, // chr7
            45200000, // chr8
            51800000, // chr9
            40300000, // chr10
            52900000, // chr11
            35400000, // chr12
            16000000, // chr13
            15600000, // chr14
            17000000, // chr15
            38200000, // chr16
            22200000, // chr17
            16100000, // chr18
            28500000, // chr19
            27100000, // chr20
            12300000, // chr21
            11800000, // chr22
            59500000, // X
            11300000, // Y
    };

    public static final Map<String, String> chrStr2Int = new HashMap<String, String>();
    public final static String doublePatternStr = "[-+]?[\\d]*\\.?[\\d]+([eE][-+]?[\\d]+)?";

    static {
        for (int i = 1; i <= 22; ++i) {
            chrStr2Int.put("" + i, "" + i);
            chrStr2Int.put("chr" + i, "" + i);
        }
        chrStr2Int.put("x", "23");
        chrStr2Int.put("X", "23");
        chrStr2Int.put("chrx", "23");
        chrStr2Int.put("chrX", "23");
        chrStr2Int.put("y", "24");
        chrStr2Int.put("Y", "24");
        chrStr2Int.put("chry", "24");
        chrStr2Int.put("chrY", "24");
    }

    @SuppressWarnings("resource")
    public static void readFile(String path) throws IOException {
        BufferedReader br = Env.getBufferedReader(path);
        for (int i = 0; i < 100; ++i) {
            System.out.println(br.readLine());
        }
    }

    public static void standardizeFormat(String inputPath, String outputPath) throws IOException {
        BufferedReader br = Env.getBufferedReader(inputPath);
        BufferedWriter bw = Env.getBufferedWriter(outputPath);

        bw.append("Chromosome\tProbeName\tPosition\tSignal\n");

        String line = br.readLine();
        while (!line.startsWith("FEATURES")) {
            line = br.readLine();
        }

        line = br.readLine();

        int probeNameIx = 8;
        int locationIx = 10;
        int logRatioIx = 14;
        Pattern p = Pattern.compile("chr([0-9XYxy]+):([0-9]+)-([0-9]+)");
        int chromosomeGroup = 1;
        int startGroup = 2;
        int endGroup = 3;
        String ratioString, ratioSign;
        String minus = "-", empty = "";

        while (line != null) {
            String[] columns = line.split("\t");
            Matcher matcher = p.matcher(columns[locationIx]);
            if (matcher.matches()) {

                bw.append(columns[probeNameIx]).append("\t");
                String chromosome = matcher.group(chromosomeGroup);
                if ("X".equalsIgnoreCase(chromosome)) {
                    chromosome = "23";
                } else if ("Y".equalsIgnoreCase(chromosome)) {
                    chromosome = "24";
                }
                bw.append(chromosome).append("\t");
                int position = (Integer.parseInt(matcher.group(startGroup)) + Integer.parseInt(matcher.group(endGroup))) / 2;
                bw.append(position + "\t");
                ratioString = prettyPrint(Double.parseDouble(columns[logRatioIx]), 2);
                bw.append(ratioString).append("\n");
            }
            line = br.readLine();
        }
        bw.close();
        br.close();
    }

    public static void toFullFormat(String inputPath, String outputPath) throws IOException {
        BufferedReader br = Env.getBufferedReader(inputPath);
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));

        bw.append("ProbeName\tChr\tStart\tEnd\tMid\tLogRatio\tLogRatioError\tPValueLogRatio\n");

        String line = br.readLine();
        while (!line.startsWith("FEATURES")) {
            line = br.readLine();
        }

        line = br.readLine();

        int probeNameIx = 8;
        int locationIx = 10;
        int logRatioIx = 14;
        int logRatioErrorIx = 15;
        int pValueLogRatioIx = 16;
        Pattern p = Pattern.compile("chr([0-9XYxy]+):([0-9]+)-([0-9]+)");
        int chromosomeGroup = 1;
        int startGroup = 2;
        int endGroup = 3;
        String ratioString, ratioSign, ratioErrorString;
        String minus = "-", empty = "";

        while (line != null) {
            String[] columns = line.split("\t");
            Matcher matcher = p.matcher(columns[locationIx]);
            if (matcher.matches()) {

                bw.append(columns[probeNameIx]).append("\t");
                String chromosome = matcher.group(chromosomeGroup);
                if ("X".equalsIgnoreCase(chromosome)) {
                    chromosome = "23";
                } else if ("Y".equalsIgnoreCase(chromosome)) {
                    chromosome = "24";
                }
                bw.append(chromosome).append("\t");

                int start = Integer.parseInt(matcher.group(startGroup));
                int end = Integer.parseInt(matcher.group(endGroup));
                int mid = (start + end) / 2;
                bw.append(start + "\t");
                bw.append(end + "\t");
                bw.append(mid + "\t");
                bw.append(prettyPrint(Double.parseDouble(columns[logRatioIx]), 5)).append("\t");
                bw.append(prettyPrint(Double.parseDouble(columns[logRatioErrorIx]), 5)).append("\t");
                bw.append(prettyPrint(Double.parseDouble(columns[pValueLogRatioIx]), 5)).append("\n");
            }
            line = br.readLine();
        }
        bw.close();
        br.close();
    }

    private static String prettyPrint(double d, int digits) {
        String signStr;
        if (d < 0) {
            signStr = "-";
            d = -d;
        } else {
            signStr = "";
        }
        int factor = (int) Math.pow(10, digits);
        int factored = (int) (d * factor + 0.5);
        String fractionStr = "" + factored % factor;
        while (fractionStr.length() < digits) {
            fractionStr = "0" + fractionStr;
        }
        return signStr + factored / factor + "." + fractionStr;
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        String sample, prefix;
        String dir = "data/cgh/";

        prefix = "US45102909_251406814090_S02_CGH-v4_91_";
        sample = "SKCO-1";

        prefix = "US45102909_251469315384_S01_CGH-v4_95_Feb07_";
        sample = "UM-UC-3";

        prefix = "US45102909_251406814400_S01_CGH-v4_91_";
        sample = "NCI-H508";


        String inputPath = dir + sample + File.separatorChar + prefix + sample + ".txt";

        //		dir = "data/CGP/";

        dir = "data/NCI-H508/";
        dir = "data/Ried/";

        //		for (double baselinePresentage = 0.01; baselinePresentage < 1.1; baselinePresentage += 0.05){
        //
        //			double[] baseline = getBaseline(dir + "NCI-H508.NCBI36_CN.CGHweb.txt", baselinePresentage, 2);
        //			System.out.println("Baseline presentage: " + baselinePresentage + ", baseline value: " + baseline[0] +", std: " + baseline[1]);
        ////			System.out.println(" (" + 2 * Math.pow(2, baseline[0] - baseline[1]) +
        ////					", " + 2 * Math.pow(2, baseline[0] + baseline[1]) + ")");
        //		}


        //		cvs2CGHweb(dir + "NCI-H508.NCBI36.genotypes.csv.gz", dir + "NCI-H508.NCBI36_SNP.CGHweb.txt");
        //		segmentation2plotFormat(dir + "NCI-H508/DNAcopy_segmentation.txt", dir + "NCI-H508/DNAcopy_segmentation.plot.txt");

        //				selectChromosome(dir + "Table_of_aCGH_smoothed_profiles.txt", dir + "NCI-H508.chr14.segmentation.txt", 14);

        //		cvs2CGHwebDir(dir);

//		inputPath = "data/Ried/";
        String outputPath = "data/Ried/segmentation/";
        double alpha = 0.001;

        String identifier = ".NCBI36.";
        //		String identifier = ".counts.";
        //		identifier = "H508";

        double stdFilterFactor = 1.5;
        double outlierFilterFactor = 0.05;
        double filteredSizeFactor = 0.3;

        double minWeightRatio = 0.2;
        double windowFraction = 0.01;

        inputPath = "data/CGP/CN/";
        outputPath = "data/CGP/CN/";
        identifier = ".csv.gz";
        File inputDir = new File(inputPath);

        String[] files = inputDir.list();
        for (int i = 0; i < files.length; ++i) {
            String file = files[i];
            if (file.contains(identifier)) {
                long start = System.currentTimeMillis();
                sample = file.substring(0, file.indexOf("."));

//				System.out.print("Standardizing " + file + "... ");
                standardizeFormat(inputPath + file, CGPFileStatistics.SEPARATOR_PATTERN_STR,
                        CGPFileStatistics.ENTRY_PATTERN_STR, "#CHR", "AFFY_ID", "POSITION",
                        "COPYNUMBER_INTENSITY", "CN_", outputPath + sample + ".CN.csv.gz");


                //				System.out.print("Analizing " + file + "... ");
                //				ChromosomeArm[][] arms = getCounts(inputPath + file, sample, minWeightRatio);

                System.out.print("Segmenting " + file + "... ");
                calcSegmentation(inputPath + file, outputPath, sample, alpha);
//
//				System.out.print("Computing counts for " + sample + "... ");
//				calcCounts(inputPath + file, outputPath + sample + ".segmentation.txt", 
//						2, windowFraction, 0.1, outputPath + sample + ".counts.txt", 
//						stdFilterFactor, outlierFilterFactor, filteredSizeFactor);
                System.out.println("Done (" + (System.currentTimeMillis() - start) / 1000 + " sec).");
            }
        }

        //		inputPath = dir + sample + ".NCBI36.genotypes.csv";
        //		printFirstLines(inputPath, 15);
        //
        ////		String oututPath = "data/cgh/wavi.txt";
        ////		toWaviFormat(inputPath, oututPath);
        ////
        //		String oututPath = dir + sample + File.separatorChar + sample + ".full.txt";
        //		toFullFormat(inputPath, oututPath);
    }

    private static void printFirstLines(String inputPath, int numOfLines) throws IOException {
        BufferedReader br = Env.getBufferedReader(inputPath);

        for (int i = 0; i < numOfLines; ++i) {
            System.out.println(br.readLine());
        }
    }

    public static void cvs2CGHwebDir(String path) throws IOException {
        File directory = new File(path);
        if (directory.isDirectory()) {
            String[] files = directory.list();
            for (String file : files) {
                int suffixStart = file.indexOf(".genotypes.csv");
                if (suffixStart > 0) {
                    cvs2CGHweb(path + file,
                            path + file.substring(0, suffixStart) + ".CN.txt");
                }
            }
        }

    }


    public static void cvs2CGHweb(String inputPath, String outputPath) throws IOException {


        BufferedReader br = Env.getBufferedReader(inputPath);
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));

        bw.append("ProbeID\tChromosome\tPosition\tLogRatio\n");
        int probeIDIx = 2;
        int chromosomeIx = 0;
        int positionIx = 1;
        int logRatioIx = 6;
        Pattern cromosomePattern = Pattern.compile("chr([\\dXYxy]+)");
        int chromosomeGroup = 1;

        String line = br.readLine();
        line = br.readLine();
        while (line != null) {
            String[] columns = line.split(",");
            if (columns[probeIDIx].startsWith("CN")) {
                Matcher matcher = cromosomePattern.matcher(columns[chromosomeIx]);
                if (matcher.matches()) {
                    bw.append(columns[probeIDIx]).append("\t");
                    String chromosome = matcher.group(chromosomeGroup);
                    if ("X".equalsIgnoreCase(chromosome)) {
                        chromosome = "23";
                    } else if ("Y".equalsIgnoreCase(chromosome)) {
                        chromosome = "24";
                    }
                    bw.append(chromosome).append("\t");
                    bw.append(columns[positionIx] + "\t");
                    bw.append(columns[logRatioIx]).append("\n");
                }
            }
            line = br.readLine();
        }
        bw.close();
        br.close();

    }

    public static void selectChromosome(String inputPath, String outputPath, int chromIx) throws IOException {
        BufferedReader br = Env.getBufferedReader(inputPath);
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));

        String chromStr = "" + chromIx;

        Pattern cromosomePattern = Pattern.compile("[^\\s]+\\s+(\\d+)");
        int chromosomeGroup = 1;
        Matcher matcher;

        bw.append(br.readLine()).append("\n");
        String line = br.readLine();
        while (line != null) {
            matcher = cromosomePattern.matcher(line);
            if (matcher.find() && chromStr.equals(matcher.group(chromosomeGroup))) {
                bw.append(line).append("\n");
            }
            line = br.readLine();
        }
        bw.close();
        br.close();
    }

    public static void segmentation2plotFormat(String inputPath, String outputPath) throws IOException {


        BufferedReader br = Env.getBufferedReader(inputPath);
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));

        bw.append("Segment\tChromosome\tPosition\tLogRatio\n");
        int segmentIx = 0;
        int chromosomeIx = 2;
        int startIx = 3;
        int endIx = 4;
        int logRatioIx = 6;

        String line = br.readLine();
        line = br.readLine();
        while (line != null) {
            String[] columns = line.split("\\s+");

            bw.append(columns[segmentIx]).append("\t");
            bw.append(columns[chromosomeIx]).append("\t");
            bw.append(columns[startIx] + "\t");
            bw.append(columns[logRatioIx]).append("\n");

            bw.append(columns[segmentIx]).append("\t");
            bw.append(columns[chromosomeIx]).append("\t");
            bw.append(columns[endIx] + "\t");
            bw.append(columns[logRatioIx]).append("\n");

            line = br.readLine();
        }
        bw.close();
        br.close();

    }

    public static double[] getBaselineSignal(String path, double windowFraction,
                                             int signalColumn, double logBase) throws IOException {
        BufferedReader br = Env.getBufferedReader(path);
        List<Double> signals = new ArrayList<Double>();

        Pattern p = Pattern.compile("([^\\s]+\\s+){" + (signalColumn - 1) + "}(" + doublePatternStr + ")");
        int logRatioGroup = 2;
        Matcher matcher;

        String line = br.readLine();
        while (line != null) {
            matcher = p.matcher(line);
            if (matcher.find()) {
                signals.add(Double.parseDouble(matcher.group(logRatioGroup)));
            }
            line = br.readLine();
        }
        br.close();

        Collections.sort(signals);

        //		//debug
        //		BufferedWriter bw1 = new BufferedWriter(new FileWriter("bbb.csv"));
        //		for (Double ratio : ratios){
        //			bw1.append("" + ratio).append("\n");
        //		}
        //		bw1.close();


        int baselineNumOfElements;
        if (windowFraction > 0 && windowFraction < 1) {
            baselineNumOfElements = (int) (1 + signals.size() * windowFraction);
        } else {
            baselineNumOfElements = signals.size();
        }


        double sum = 0, sumOfSquares = 0;

        double average, variance, signal;

        double baselineValue = 0;
        double baselineVariance = Double.MAX_VALUE;

        ////		debug
        //		BufferedWriter bw2 = new BufferedWriter(new FileWriter("ccc.csv"));

        //		double[] sumVec = new double[ratios.size()+1];
        //		double[] sumOfSquaresVec = new double[ratios.size()+1];
        //
        //		for (int i=0; i<ratios.size(); ++i){
        //			ratio = ratios.get(i);
        //			sumVec[i+1] = sumVec[i] + ratio;
        //			sumOfSquaresVec[i+1] = sumOfSquaresVec[i] + ratio * ratio;
        //		}
        //
        //		double minWindowFraction = 0.01, maxWindowFraction = 0.8;
        //		int segments = 10;
        //		for (int j = 0; j < 6; ++j){
        //			int bestSegmentStart = bestStdWindowLength(sumVec, sumOfSquaresVec, minWindowFraction, maxWindowFraction, segments);
        //			double step = (maxWindowFraction - minWindowFraction) / segments;
        //			maxWindowFraction = Math.min(1, minWindowFraction + step*(bestSegmentStart+1));
        //			minWindowFraction = Math.max(0.01, minWindowFraction + step*(bestSegmentStart-1));
        //		}


        int i = 0;
        for (; i < baselineNumOfElements - 1; ++i) {
            signal = signals.get(i);
            sum += signal;
            sumOfSquares += signal * signal;
        }

        for (; i < signals.size(); ++i) {
            signal = signals.get(i);
            sum += signal;
            sumOfSquares += signal * signal;

            average = sum / baselineNumOfElements;

            variance = sumOfSquares + average * (baselineNumOfElements * average - 2 * sum);

            if (baselineVariance > variance) {
                baselineVariance = variance;
                baselineValue = average;
                //				bw2.append(baselineValue + "\t,\t" + (i-baselineNumOfElements/2)).append("\n");
            }
            signal = signals.get(i - baselineNumOfElements + 1);
            sum -= signal;
            sumOfSquares -= signal * signal;
        }

        //		bw2.close();


        double std = Math.sqrt(baselineVariance / (baselineNumOfElements - 1));

        return new double[]{baselineValue, std};

    }

    private static int bestStdWindowLength(double[] sumVec,
                                           double[] sumOfSquaresVec, double minWidowFraction,
                                           double maxWidowFraction, int segments) {
        double step = (maxWidowFraction - minWidowFraction) / segments;
        double bestVar = Double.POSITIVE_INFINITY;
        int bestWindowIx = -1;
        int bestI = -1;

        int windowIx = 0;
        for (double widowFraction = minWidowFraction; widowFraction <= maxWidowFraction; widowFraction += step, ++windowIx) {
            int widowSize = (int) (widowFraction * (sumVec.length - 1));
            for (int i = widowSize; i < sumVec.length; ++i) {
                double average = (sumVec[i] - sumVec[i - widowSize]) / widowSize;
                double currVar = (sumOfSquaresVec[i] - sumOfSquaresVec[i - widowSize] - widowSize * average * average) / (widowSize - 1);
                if (bestVar > currVar) {
                    bestVar = currVar;
                    bestWindowIx = windowIx;
                    bestI = i - widowSize / 2;
                }
            }
        }

        return bestWindowIx;
    }

    public static void calcSegmentation(String inputPath, String outputPath, String id, double alpha) throws InterruptedException {

        try {
            // Execute command
//			"Chromosome" + CSV_SEPARATOR + "ProbeName" + CSV_SEPARATOR + "Position" + CSV_SEPARATOR + "Signal
            String command = "Rscript R/segmentCGH.R " + inputPath + " " + outputPath + " Signal Chromosome Position " + id + " " + alpha; // + " > " + outputPath;
            Process child = Runtime.getRuntime().exec(command);
            child.waitFor();
        } catch (IOException e) {
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
    }

    public static void calcCounts(String cghFile, String segmentationFile,
                                  double logBase, double windowFraction, double trimProp, String countsFile,
                                  double stdFilterFactor, double outlierFilterFactor,
                                  double filteredSizeFactor) throws IOException {
        int signalColumn = 4;


        double[] baseline = getBaselineSignal(cghFile, windowFraction, signalColumn, logBase);

        //		double factor = 2/baseline[0];
        double factor = Math.log(logBase); // transforms between logBase and the natural logarithm base e.
        double addativeFactor = Math.log(2) / factor - baseline[0];

        //		double medianError1 = 0, logMedianError1 = 0, medianError2 = 0, logMedianError2 = 0;

        BufferedReader segmentationReader = Env.getBufferedReader(segmentationFile);
        BufferedReader cghReader = Env.getBufferedReader(cghFile);
        List<Double> signals = new ArrayList<Double>();
        //		List<Double> countWeights = new ArrayList<Double>();

        BufferedWriter countsWriter = new BufferedWriter(new FileWriter(countsFile));
        //		countsWriter.append("Segment\tChromosome\tStart\tEnd\tMarkers\tMinCount\tMaxCount\tBestCount\tProbabilities\n");
        countsWriter.append("Segment").append(CSV_SEPARATOR);
        countsWriter.append("Chromosome").append(CSV_SEPARATOR);
        countsWriter.append("Start").append(CSV_SEPARATOR);
        countsWriter.append("End").append(CSV_SEPARATOR);
        countsWriter.append("Signals").append(CSV_SEPARATOR);
        countsWriter.append("NonFilteredSignals").append(CSV_SEPARATOR);
        countsWriter.append("SumOfSignals").append(CSV_SEPARATOR);
        countsWriter.append("SumOfSquareSignals").append(CSV_SEPARATOR);
        countsWriter.append("SumOfLogSignals").append(CSV_SEPARATOR);
        countsWriter.append("SumOfLogSquareSignals").append("\n");


        Pattern cghPattern = Pattern.compile("[^\\s]+\\s+([\\d]+)\\s+([\\d]+)\\s+(" + doublePatternStr + ")");
        int cghChrGroup = 1;
        int cghPosGroup = 2;
        int cghRatioGroup = 3;

        Pattern segPattern = Pattern.compile("([\\d]+)\\s+[^\\s]+\\s+([\\d]+)\\s+([\\d]+)\\s+([\\d]+)\\s+([\\d]+)\\s+");
        int segmentGroup = 1;
        int segChrGroup = 2;
        int startGroup = 3;
        int endGroup = 4;
        //		int markersGroup = 5;

        Matcher matcher;

        String cghChr = null, segChr;
        int cghPos = 0, segStart, segEnd;
        double signal = 0;

        String cghLine = cghReader.readLine();
        while (cghLine != null) {
            matcher = cghPattern.matcher(cghLine);
            if (matcher.find()) {
                cghChr = matcher.group(cghChrGroup);
                cghPos = Integer.parseInt(matcher.group(cghPosGroup));
                //				signal = Double.parseDouble(matcher.group(cghRatioGroup)); // factor * Math.pow(logBase, Double.parseDouble(matcher.group(cghRatioGroup)));
                signal = factor * (addativeFactor + Double.parseDouble(matcher.group(cghRatioGroup)));
                break;
            }
            cghLine = cghReader.readLine();
        }

        if (cghLine != null) {
            String segLine = segmentationReader.readLine();
            while (segLine != null) {
                matcher = segPattern.matcher(segLine);
                if (matcher.find()) {
                    String segment = matcher.group(segmentGroup);
                    segChr = matcher.group(segChrGroup);
                    segStart = Integer.parseInt(matcher.group(startGroup));
                    segEnd = Integer.parseInt(matcher.group(endGroup));
                    signals.clear();

                    while (cghLine != null && segChr.equals(cghChr) && segStart <= cghPos && segEnd >= cghPos) {
                        signals.add(signal);
                        cghLine = cghReader.readLine();
                        if (cghLine != null) {
                            matcher = cghPattern.matcher(cghLine);
                            if (matcher.find()) {
                                cghChr = matcher.group(cghChrGroup);
                                cghPos = Integer.parseInt(matcher.group(cghPosGroup));
                                //								signal = Double.parseDouble(matcher.group(cghRatioGroup)); // factor * Math.pow(logBase, Double.parseDouble(matcher.group(cghRatioGroup)));
                                signal = addativeFactor + factor * Double.parseDouble(matcher.group(cghRatioGroup));
                            }
                        }
                    }

                    if (!signals.isEmpty()) {
                        Collections.sort(signals);
                        //						int minCount = (int) ((double) signals.get(0));
                        //						int bestCount = -1;
                        //						int currCount = minCount;
                        //						double bestCountWeight = 0;
                        //						double currCountWeight = 0;
                        //						double nextCountWeight = 0;
                        //						countWeights.clear();
                        int numOfSignals = signals.size();
                        //						signals.add(signals.get(numOfSignals-1) + 1.1);

                        double sum = 0, sumOfSquares = 0, sumOfLogs = 0, sumOfSquareLogs = 0;


                        for (int i = 0; i < numOfSignals; ++i) {
                            double currSignal = signals.get(i);
                            double currSignalExp = Math.exp(currSignal);

                            sumOfLogs += currSignal;
                            sumOfSquareLogs += currSignal * currSignal;
                            sum += currSignalExp;
                            sumOfSquares += currSignalExp * currSignalExp;


                            //							while (currSignal - currCount > 1){
                            //								countWeights.add(currCountWeight);
                            //								if (bestCountWeight < currCountWeight){
                            //									bestCountWeight = currCountWeight;
                            //									bestCount = currCount;
                            //								}
                            //								if ((numOfSignals-i+nextCountWeight)/bestCountWeight < trimProp && !countWeights.isEmpty()){
                            //									i = numOfSignals;
                            //									break;
                            //								}
                            //								++currCount;
                            //								currCountWeight = nextCountWeight;
                            //								nextCountWeight = 0;
                            //							}
                            //							currCountWeight += currCount + 1 - currSignal;
                            //							nextCountWeight += currSignal - currCount;
                        }


                        //						int trimStart = 0, trimEnd = countWeights.size()-1;
                        //						for (; countWeights.get(trimStart)/bestCountWeight < trimProp; ++trimStart);
                        //						for (; countWeights.get(trimEnd)/bestCountWeight < trimProp; --trimEnd);


                        double mean = sum / numOfSignals;
                        double var = (sumOfSquares - mean * mean * numOfSignals) / (numOfSignals - 1);
                        double std = Math.sqrt(var);
                        double logMean = sumOfLogs / numOfSignals;
                        double logVar = (sumOfSquareLogs + logMean * logMean * numOfSignals) / (numOfSignals - 1);
                        double logStd = Math.sqrt(logVar);

                        //						int midIx = numOfSignals/2;
                        //						double med1 = signals.get(midIx);
                        //						double med2 = signals.get(numOfSignals-midIx-1);
                        //						double median = (Math.exp(med1) + Math.exp(med2))/2;
                        //						double expectedLogNormalMedian = Math.exp(logMean);
                        //						double logMedianErr1 = Math.pow(1-expectedLogNormalMedian/median, 2);
                        //						double medianErr1 = Math.pow(1-mean/median, 2);
                        //						double logMedianErr2 = (median - expectedLogNormalMedian)/ median;
                        //						double medianErr2 = (median - mean)/ median;
                        //
                        //						medianError1 += medianErr1;
                        //						medianError2 += medianErr2;
                        //						logMedianError1 += logMedianErr1;
                        //						logMedianError2 += logMedianErr2;

                        if (logStd > stdFilterFactor * logMean) {
                            System.err.println("Filtering segment: " + segLine);
                            System.err.println("Segment count mean: " + logMean + ", segment count std: " + logStd + ", ratio: " + (logStd / logMean));
                            continue;
                        }

                        NormalDistribution nd = new NormalDistribution(logMean, logStd);

                        int leftMarginIx;
                        double outlyerDev = logMean - nd.inverseCumulativeProbability(outlierFilterFactor);
                        //						if (outlyerDev < logMean){
                        leftMarginIx = Collections.binarySearch(signals, logMean - outlyerDev);
                        if (leftMarginIx < 0) {
                            leftMarginIx = -leftMarginIx - 1;
                        }
                        //						}
                        else {
                            leftMarginIx = 0;
                        }


                        double rightMarginIx = Collections.binarySearch(signals, Math.log(mean + outlyerDev));
                        if (rightMarginIx < 0) {
                            rightMarginIx = -rightMarginIx - 1;
                        }

                        if (1 - (rightMarginIx - leftMarginIx) / signals.size() > filteredSizeFactor) {
                            System.err.println("Filtering segment: " + segLine);
                            System.err.println("Proprtion of filtered signals: " + (1 - (rightMarginIx - leftMarginIx) / signals.size()));
                            continue;
                        }

                        for (int i = 0; i < leftMarginIx; ++i) {
                            double currSignal = signals.get(i);
                            double currSignalExp = Math.exp(currSignal);

                            sumOfLogs -= currSignal;
                            sumOfSquareLogs -= currSignal * currSignal;
                            sum -= currSignalExp;
                            sumOfSquares -= currSignalExp * currSignalExp;
                        }

                        for (int i = signals.size() - 1; i >= rightMarginIx; --i) {
                            double currSignal = signals.get(i);
                            double currSignalExp = Math.exp(currSignal);

                            sumOfLogs -= currSignal;
                            sumOfSquareLogs -= currSignal * currSignal;
                            sum -= currSignalExp;
                            sumOfSquares -= currSignalExp * currSignalExp;
                        }

                        numOfSignals = (int) (rightMarginIx - leftMarginIx);

                        countsWriter.append(segment).append(CSV_SEPARATOR);
                        countsWriter.append(segChr).append(CSV_SEPARATOR);
                        countsWriter.append("" + segStart).append(CSV_SEPARATOR);
                        countsWriter.append("" + segEnd).append(CSV_SEPARATOR);
                        countsWriter.append("" + signals.size()).append(CSV_SEPARATOR);
                        countsWriter.append("" + numOfSignals).append(CSV_SEPARATOR);
                        countsWriter.append("" + sum).append(CSV_SEPARATOR);
                        countsWriter.append("" + sumOfSquares).append(CSV_SEPARATOR);
                        countsWriter.append("" + sumOfLogs).append(CSV_SEPARATOR);
                        countsWriter.append("" + sumOfSquareLogs).append("\n");

                        //
                        //						mean = sum/numOfSignals;
                        //						var = (sumOfSquares - mean*mean*numOfSignals) / (numOfSignals-1);
                        //						std = Math.sqrt(var);
                        //						logMean = sumOfLogs/numOfSignals;
                        //						logVar = (sumOfSquareLogs + logMean*logMean*numOfSignals) / (numOfSignals-1);
                        //						logStd = Math.sqrt(logVar);
                        //
                        //
                        //						double logNormalMean = Math.exp(logMean+logVar/2);
                        //						double logNormalMed = Math.exp(logMean);
                        //						double logNormalMode = Math.exp(logMean-logVar);
                        //
                        //						countsWriter.append(segment + "\t" + segChr + "\t" +
                        //								segStart + "\t" + segEnd + "\t" + numOfSignals + "\t");
                        //						double error = 0, logError = 0;
                        //						for (int i = 0; i < numOfSignals; ++i){
                        //							double currSignal = signals.get(i);
                        //							double currSignalExp = Math.exp(currSignal);
                        //							error += (currSignalExp - mean) * (currSignalExp - mean);
                        //							logError += (currSignalExp - logNormalMean) * (currSignalExp - logNormalMean);
                        //						}
                        //
                        //						countsWriter.append(mean + "\t" + var + "\t");
                        //						countsWriter.append(logMean + "\t" + logVar + "\t");
                        //						countsWriter.append(logNormalMean + "\t" + logNormalMed
                        //								+ "\t" + logNormalMode + "\t" + (Math.exp(logVar)-1) * Math.exp(2*logMean+logVar));
                        //
                        //
                        ////								+ (minCount+trimStart) +"\t" + (minCount+trimEnd)
                        ////								+"\t" + bestCount + "\t");
                        ////						for (int i=trimStart; i<=trimEnd; ++i){
                        ////							countsWriter.append(countWeights.get(i)/numOfSignals + " ");
                        ////						}
                        //						countsWriter.append("\n");
                    }
                }
                segLine = segmentationReader.readLine();
            }
        }
        segmentationReader.close();
        cghReader.close();
        countsWriter.close();

        //		System.err.println(medianError1 + ", " + medianError2 + ", " + logMedianError1 + ", " + logMedianError2);
    }

    public static List<ChromosomeArm> getAllChromosomeArms(String path, String fileIdentifier, double minWeightRatio, double weightFactor) throws IOException {
        List<ChromosomeArm> arms = new ArrayList<ChromosomeArm>();
        String sample;
        File inputDir = new File(path);
        for (String file : inputDir.list()) {
            if (file.contains(fileIdentifier)) {
                //				long start = System.currentTimeMillis();
                sample = file.substring(0, file.indexOf("."));
                //				System.out.print("Analizing " + file + "... ");
                ChromosomeArm[][] sampleArms = getCounts(path + file, sample, minWeightRatio, weightFactor);

                //				System.out.print("Segmenting " + file + "... ");
                //				//				calcSegmentation(inputPath + file, outputPath, sample, alpha);
                //				System.out.print("Computing counts for " + sample + "... ");
                //				calcCounts(inputPath + file, outputPath + sample + ".segmentation.txt",
                //						2, 0.3, 0.1, outputPath + sample + ".counts.txt",
                //						stdFilterFactor, outlierFilterFactor, filteredSizeFactor);
                //				System.out.println("Done (" + (System.currentTimeMillis()-start)/1000 + " sec).");

                for (int i = 0; i < sampleArms.length; ++i) {
                    for (int j = 0; j < sampleArms[i].length; ++j) {
                        if (sampleArms[i][j] != null) {
                            arms.add(sampleArms[i][j]);
                        }
                    }
                }
            }
        }

        return arms;
    }

    public static ChromosomeArm[][] getCounts(String countsFile, String sampleName, double minWeightRation, double weightFactor) throws IOException {
        ChromosomeArm[][] arms = new ChromosomeArm[25][2];

        BufferedReader br = Env.getBufferedReader(countsFile);
        List<Integer> minCounts = new ArrayList<Integer>();
        List<List<Double>> probablities = new ArrayList<List<Double>>();
        List<Integer> start = new ArrayList<Integer>();
        List<Integer> end = new ArrayList<Integer>();

        String gIntP = "([\\d]+)\\t";
        String intP = "[\\d]+\\t";
        String gDoubleP = "(" + doublePatternStr + ")\\t";

        Pattern p = Pattern.compile(intP + gIntP + gIntP + gIntP + intP + gIntP + gDoubleP + gDoubleP + gDoubleP + gDoubleP);
        int chrGroup = 1;
        int startGroup = 2;
        int endGroup = 3;
        int signalsGroup = 4;
        int sumGroup = 5;
        int sumOfSquaresGroup = 7;
        int sumOfLogsGroup = 9;
        int sumOfLogSquaresGroup = 11;

        Matcher matcher;

        String chrStr = null;
        int chrIx, segStart, segEnd, signals, currCount, countMedian;
        double sum, sumOfSquares, average, maxWeight, currLog, logVar, var, std;

        int currChromosome = -1;
        boolean isP = true;
        ChromosomeArm currArm = null;

        String line = br.readLine();
        while (line != null) {
            line += "\t";
            matcher = p.matcher(line);
            if (matcher.find()) {
                chrStr = matcher.group(chrGroup);
                chrIx = Integer.parseInt(chrStr);
                segStart = Integer.parseInt(matcher.group(startGroup));
                segEnd = Integer.parseInt(matcher.group(endGroup));
                signals = Integer.parseInt(matcher.group(signalsGroup));
                sum = Double.parseDouble(matcher.group(sumOfLogsGroup));
                sumOfSquares = Double.parseDouble(matcher.group(sumOfLogSquaresGroup));
                average = sum / signals;
                countMedian = (int) Math.exp(average);
                double lowerLog = Math.log(countMedian);
                double upperLog = Math.log(countMedian + 1);
                if (average - lowerLog < upperLog - average) {
                    currLog = lowerLog;
                } else {
                    currLog = upperLog;
                    ++countMedian;
                }

                logVar = (sumOfSquares + currLog * (signals * currLog - 2 * sum)) / (signals - 1);
                var = (Math.exp(logVar) - 1) * Math.exp(2 * currLog + logVar);
                std = Math.sqrt(var);
                maxWeight = Math.exp(-std * weightFactor);
                if (maxWeight == 0) {
                    // segment is ignored due to high variance
                    line = br.readLine();
                    continue;
                }

                double currWeight = maxWeight;
                List<Double> currProbs = new ArrayList<Double>();

                for (int i = countMedian - 1; i >= -1 && currWeight / maxWeight >= minWeightRation; --i) {
                    currProbs.add(currWeight);
                    currLog = Math.log(i);
                    logVar = (sumOfSquares + currLog * (signals * currLog - 2 * sum)) / (signals - 1);
                    var = (Math.exp(logVar) - 1) * Math.exp(2 * currLog + logVar);
                    std = Math.sqrt(var);
                    currWeight = Math.exp(-std * weightFactor);
                }
                int minCount = countMedian - currProbs.size() + 1;
                Collections.reverse(currProbs);
                currLog = Math.log(countMedian + 1);
                logVar = (sumOfSquares + currLog * (signals * currLog - 2 * sum)) / (signals - 1);
                var = (Math.exp(logVar) - 1) * Math.exp(2 * currLog + logVar);
                std = Math.sqrt(var);
                currWeight = Math.exp(-std * weightFactor);
                for (int i = countMedian + 2; currWeight / maxWeight >= minWeightRation; ++i) {
                    currProbs.add(currWeight);
                    currLog = Math.log(i);
                    logVar = (sumOfSquares + currLog * (signals * currLog - 2 * sum)) / (signals - 1);
                    var = (Math.exp(logVar) - 1) * Math.exp(2 * currLog + logVar);
                    std = Math.sqrt(var);
                    currWeight = Math.exp(-std * weightFactor);
                }

                if (currChromosome != chrIx || (isP && centromerPositions[chrIx] < segEnd)) {
                    if (currChromosome == chrIx) {
                        if (centromerPositions[chrIx] > segStart) {
                            // current segment spans the centromer
                            minCounts.add(minCount);
                            probablities.add(currProbs);
                            start.add(segStart);
                            end.add(centromerPositions[chrIx]);
                        }

                        // finalizing the p-arm of this chromosome:
                        currArm = new ChromosomeArm(sampleName, currChromosome, true, null);
                        currArm.setCounts(minCounts, probablities, start, end);
                        arms[chrIx][0] = currArm;
                        isP = false;
                    } else {
                        // finalizing the arm of the previous chromosome:
                        if (currChromosome >= 0) {
                            currArm = new ChromosomeArm(sampleName, currChromosome, isP, null);
                            currArm.setCounts(minCounts, probablities, start, end);
                            if (isP) {
                                arms[currChromosome][0] = currArm;
                            } else {
                                arms[currChromosome][1] = currArm;
                            }
                        }
                        currChromosome = chrIx;
                        isP = centromerPositions[chrIx] > segStart;
                    }

                    // clearing the lists for the next arm
                    minCounts.clear();
                    probablities.clear();
                    start.clear();
                    end.clear();
                }
                minCounts.add(minCount);
                probablities.add(currProbs);
                if (!isP && centromerPositions[chrIx] > segStart) {
                    start.add(centromerPositions[chrIx]);
                } else {
                    start.add(segStart);
                }
                end.add(segEnd);
            }
            line = br.readLine();
        }

        currArm = new ChromosomeArm(sampleName, currChromosome, isP, null);
        currArm.setCounts(minCounts, probablities, start, end);
        if (isP) {
            arms[currChromosome][0] = currArm;
        } else {
            arms[currChromosome][1] = currArm;
        }

        return arms;
    }


    public static void standardizeFormat(String inputPath, String sepPtrnStr, String entPtrnStr,
                                         String chromosomeTitle, String probeNameTitle, String locationTitle,
                                         String signalTitle, String validProbeIdentifier, String outputPath) throws IOException {

        BufferedReader br = Env.getBufferedReader(inputPath);
        BufferedWriter bw = Env.getBufferedWriter(outputPath);

        bw.append("Chromosome" + CSV_SEPARATOR + "ProbeName" + CSV_SEPARATOR + "Position" + CSV_SEPARATOR + "Signal\n");

        String line = br.readLine(); // titles line
        Map<String, Integer> titlesMap = Env.columnTitles2Ixs(line, sepPtrnStr);

        int chromosomeIx = titlesMap.get(chromosomeTitle);
        int locationIx = titlesMap.get(locationTitle);
        int probeNameIx = titlesMap.get(probeNameTitle);
        int signalIx = titlesMap.get(signalTitle);

        Pattern p = Env.getPattern(Env.max(chromosomeIx, locationIx, probeNameIx, signalIx), sepPtrnStr, entPtrnStr);

        line = br.readLine();
        while (line != null) {
            Matcher matcher = p.matcher(line);
            if (matcher.find()) {
                if (matcher.group(probeNameIx).contains(validProbeIdentifier)) {
                    bw.append(chrStr2Int.get(matcher.group(chromosomeIx))).append(CSV_SEPARATOR);
                    bw.append(matcher.group(probeNameIx)).append(CSV_SEPARATOR);
                    bw.append(matcher.group(locationIx)).append(CSV_SEPARATOR);
                    bw.append(matcher.group(signalIx)).append("\n");
                }
            }
            line = br.readLine();
        }
    }
}
