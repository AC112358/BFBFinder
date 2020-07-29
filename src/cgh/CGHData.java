package cgh;

import bfbf.BFBCalculator;
import bfbf.BFB_Algo;
import bfbf.Env;
import bfbf.Solution1;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class CGHData implements Serializable {

    public static final int P = 0;
    public static final int Q = 1;
    public static final double epsilone = Math.pow(10, -10);
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
    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private static final int MAX_REFINEMENTS = 10;
    private static final int MAX_REFINEMENT_WINDOW = 20;
    private static final Map<String, Integer> chrStr2Int = new HashMap<>();

    static {
        for (int i = 1; i <= 22; ++i) {
            chrStr2Int.put("" + i, i);
            chrStr2Int.put("chr" + i, i);
        }
        chrStr2Int.put("x", 23);
        chrStr2Int.put("X", 23);
        chrStr2Int.put("chrx", 23);
        chrStr2Int.put("chrX", 23);
        chrStr2Int.put("y", 24);
        chrStr2Int.put("Y", 24);
        chrStr2Int.put("chry", 24);
        chrStr2Int.put("chrY", 24);
    }


    protected String sampleName;
    protected int[][][] positions;

    /**
     * lnSignals holds the lan of the measured probe intensity. It is assume
     * that the value given in the file is of the form y = log_b(x/2), where b
     * is the log base and x/2 is the measured red/green ratio. Since we assume
     * the reference genome has two copies of the segment probe, the kept value
     * is y * ln(b) + ln(2) = ln(x).
     */
    protected float[][][] lnSignals;
    protected float[][][] lnSignalSums;
    protected float[][][] lnSignalSquareSums;
    protected TIntList[][] segStarts;
    protected TIntList[][] segEnds;

    protected double shiftFactor;

    public CGHData(String cgpFilePath) throws IOException {
        this(cgpFilePath, "\\s*,\\s*", "#CHR", "AFFY_ID", "POSITION",
                "COPYNUMBER_INTENSITY", "CN_", "SEGMENT_COPYNUMBER", 2);
    }

    public CGHData(String path, String sepPtrnStr,
                   String chromosomeTitle, String probeNameTitle, String locationTitle,
                   String signalTitle, String validProbeIdentifier, String segmentTitle, double logBase) throws IOException {

        BufferedReader br = Env.getBufferedReader(path);
        String line = br.readLine(); // titles line
        Map<String, Integer> titlesMap = Env.columnTitles2Ixs(line, sepPtrnStr);

        boolean readSegments = segmentTitle != null;
        int chromosomeIx = titlesMap.get(chromosomeTitle) - 1;
        int positionIx = titlesMap.get(locationTitle) - 1;
        int probeNameIx = titlesMap.get(probeNameTitle) - 1;
        int signalIx = titlesMap.get(signalTitle) - 1;
        int segmentIx = -1;
        if (readSegments) {
            segmentIx = titlesMap.get(segmentTitle) - 1;
        }
        int maxIx = Env.max(chromosomeIx, positionIx, probeNameIx, signalIx, segmentIx);

        double logFactor = Math.log(logBase); // for transforming log base to e
        double logShift = Math.log(2); // for transforming ln(x/2) to ln(x)

        positions = new int[25][2][];
        lnSignals = new float[25][2][];
        lnSignalSums = new float[25][2][];
        lnSignalSquareSums = new float[25][2][];
        segStarts = new TIntList[25][2];
        segEnds = new TIntList[25][2];

        TIntList pos = new TIntArrayList();
        TFloatList signals = new TFloatArrayList();
        TIntList starts = new TIntArrayList();
        TIntList ends = new TIntArrayList();

        int chrIx, position;
        int currChromosome = 0;
        String currSegmentStr = "noSeg";
        int arm = P;

        line = br.readLine();
        while (line != null) {
            String[] splitted = line.split(sepPtrnStr);
            if (splitted.length >= maxIx && splitted[probeNameIx].contains(validProbeIdentifier)) {
                chrIx = chrStr2Int.get(splitted[chromosomeIx]);
                position = Integer.parseInt(splitted[positionIx]);

                if (chrIx != currChromosome || (arm == P && centromerPositions[chrIx] < position)) {
                    // finalizing the current accumulated arm
                    positions[currChromosome][arm] = pos.toArray();
                    lnSignals[currChromosome][arm] = signals.toArray();
                    lnSignalSums[currChromosome][arm] = prefixSum(lnSignals[currChromosome][arm]);
                    lnSignalSquareSums[currChromosome][arm] = prefixSumOfSquares(lnSignals[currChromosome][arm]);

                    if (readSegments) {
                        ends.add(pos.size());
                        segStarts[currChromosome][arm] = starts;
                        segEnds[currChromosome][arm] = ends;
                        starts = new TIntArrayList();
                        starts.add(0);
                        ends = new TIntArrayList();
                        currSegmentStr = splitted[segmentIx];
                    }

                    pos.clear();
                    signals.clear();
                    currChromosome = chrIx;
                    if (centromerPositions[chrIx] < position) {
                        arm = Q;
                    } else {
                        arm = P;
                    }
                } else if (readSegments && !pos.isEmpty() && !currSegmentStr.equals(splitted[segmentIx])) {
                    // adding a new segment
                    ends.add(pos.size());
                    starts.add(pos.size());
                    currSegmentStr = splitted[segmentIx];
                }

                pos.add(position);
                signals.add((float) (logFactor * Double.parseDouble(splitted[signalIx]) + logShift));
            }
            line = br.readLine();
        }

        shiftFactor = 0;
        refineShiftFactor();
        //		while(refineSegments()){
        //			refineShiftFactor();
        //		}
    }

    public static void main(String[] args) throws IOException, ClassNotFoundException {
        String path = "data/CGP/";
        File dir = new File(path);
        long start = System.currentTimeMillis();
        List<Solution1> solutions = new ArrayList<>();
        TIntList minCounts = new TIntArrayList();
        List<TDoubleList> countProbs = new ArrayList<>();
        double maxError = 0.2;
        int beta = 100;

        int i = 0;
        for (String file : dir.list()) {
            if (file.contains("obj.gz")) {
                ++i;
//				if (file.contains("NCBI36.genotypes.csv.gz")){
//				String inPath = file.substring(0, file.indexOf('.')) + ".obj";
                System.out.print(i + ") Analyzing " + file + ", ");
                long currStart = System.currentTimeMillis();
//				CGHData data = new CGHData(path + file);
//				data.writeToFile(path + inPath);
//				long end = System.currentTimeMillis();
//				System.out.println("time: " + (end-currStart)/1000 + " (" + (end-start)/1000 + ")");

                CGHData data = readFromFile(path + file);
//				System.out.println("Time: " + (System.currentTimeMillis()-start)/1000);
//				System.out.println("Shift factor: " + data.shiftFactor + ", shift factor exponent: " + Math.exp(data.shiftFactor));
                //		data.printCountHistogram();

                for (int chromosome = 1; chromosome <= 24; ++chromosome) {
                    boolean isP = true;
                    for (int j = 0; j < 2; ++j) {
                        ChromosomeArm arm = data.getChromosomeArm(chromosome, isP);
                        if (arm == null) {
                            continue;
                        }
                        arm.getCountWeights(minCounts, countProbs, 1 - maxError, beta);
//						System.out.println(arm);

                        double[][] countProbsArr = new double[countProbs.size()][];
                        for (int l = 0; l < countProbs.size(); ++l) {
                            countProbsArr[l] = countProbs.get(l).toArray();
                        }

                        int[] minCountArr = minCounts.toArray();
                        BFB_Algo.fixPloidy(minCountArr, 2);
                        Solution1 res = BFBCalculator.longestBFBSubstring(minCountArr, countProbsArr, maxError);
                        if (res.getEffectiveLength() >= 6) {
                            solutions.add(res);
//							System.out.println(res);
//							System.exit(0);
                        }
//						System.out.println(res);
//						System.out.println();
                        isP = false;
                    }
                }
                long end = System.currentTimeMillis();
                System.out.println("time: " + (end - currStart) / 1000 + " (" + (end - start) / 1000 + "), found bfb vectors: " + solutions.size());
            }
        }

        for (Solution1 solution : solutions) {
            System.out.println(solution);
        }


//		String path = "data/CGP/NCI-H508.NCBI36.genotypes.csv.gz";
//		CGHData data = new CGHData(path);
//		System.out.println("Time: " + (System.currentTimeMillis()-start)/1000);
//		System.out.println("Shift factor: " + data.shiftFactor + ", shift factor exponent: " + Math.exp(data.shiftFactor));
//		//		data.printCountHistogram();
//		TIntList minCounts = new TIntArrayList();
//		List<TDoubleList> countProbs = new ArrayList<TDoubleList>();
//		double maxError = 0.3;
//		int beta = 3;
//
//		for (int chromosome = 1; chromosome <= 24; ++chromosome){
//			boolean isP = true;
//			for (int j=0; j<2; ++j){
//				ChromosomeArm arm = data.getChromosomeArm(chromosome, isP);
//				if (arm == null){
//					continue;
//				}
//				arm.getCountWeights(minCounts, countProbs, 1 - maxError, beta);
//				System.out.println(arm);
//
//				double[][] countProbsArr = new double[countProbs.size()][];
//				for (int i=0; i<countProbs.size(); ++i){
//					countProbsArr[i] = countProbs.get(i).toArray();
//				}
//				Solution1 res = BFBCalculator.longestBFBSubstring(minCounts.toArray(), countProbsArr, maxError);
//				System.out.println(res);
//				System.out.println();
//				isP = false;
//			}
//		}
//
//
        System.out.println("Done, " + solutions.size() + " BFB vectors");
    }

    @SuppressWarnings("resource")
    public static CGHData readFromFile(String path) throws IOException, ClassNotFoundException {
        return (CGHData) new ObjectInputStream(new GZIPInputStream(new FileInputStream(path))).readObject();
    }

    // TODO: accommodate corrections
    private static double[] getVars(TIntList segStarts, TIntList segEnds, float[] sums, float[] sumsOfSquares) {
        double[] sds = new double[segStarts.size()];
        for (int i = 0; i < segStarts.size(); ++i) {
            int n = segEnds.get(i) - segStarts.get(i);
            float mean = (sums[segEnds.get(i)] - sums[segStarts.get(i)]) / n;
            sds[i] = (sumsOfSquares[segEnds.get(i)] - sumsOfSquares[segStarts.get(i)]) / n - mean * mean;
        }
        return sds;
    }

    private static double correct(double mean) {
        return Math.log(Math.round(Math.exp(mean)));
    }

    private static int[] getLengths(TIntList segStarts, TIntList segEnds) {
        int[] lengths = new int[segStarts.size()];
        for (int i = 0; i < segStarts.size(); ++i) {
            lengths[i] = segEnds.get(i) - segStarts.get(i);
        }
        return lengths;
    }

    private static int[] getValues(TIntList segStarts, int[] positions, int shift) {
        int[] values = new int[segStarts.size()];
        for (int i = 0; i < segStarts.size(); ++i) {
            values[i] = positions[segStarts.get(i) + shift];
        }
        return values;
    }

    private static final int getCount(double mean) {
        return (int) Math.round(Math.exp(mean));
    }

    private void refineShiftFactor() {
        double delta = localIntegerCountFactor();
        for (int i = 0; i < MAX_REFINEMENTS && (delta > epsilone || delta < -epsilone); ++i) {
            shiftFactor += delta;
            delta = localIntegerCountFactor();
        }
    }

    private float[] prefixSum(float[] signals) {
        float[] sums = new float[signals.length + 1];
        for (int i = 0; i < signals.length; ++i) {
            sums[i + 1] = sums[i] + signals[i];
        }
        return sums;
    }

    private float[] prefixSumOfSquares(float[] signals) {
        float[] sums = new float[signals.length + 1];
        for (int i = 0; i < signals.length; ++i) {
            sums[i + 1] = sums[i] + signals[i] * signals[i];
        }
        return sums;
    }

    private float[] prefixSumOfExponents(float[] signals) {
        float[] sums = new float[signals.length + 1];
        for (int i = 0; i < signals.length; ++i) {
            sums[i + 1] = (float) (sums[i] + Math.exp(signals[i]));
        }
        return sums;
    }

    public ChromosomeArm getChromosomeArm(int chromosomeIx, boolean isP) {
        int side;
        if (isP) {
            side = P;
        } else {
            side = Q;
        }

        if (segStarts[chromosomeIx][side] == null) {
            return null;
        }

        ChromosomeArm arm = new ChromosomeArm(sampleName, chromosomeIx, isP,
                getValues(segStarts[chromosomeIx][side], positions[chromosomeIx][side], 0),
                getValues(segEnds[chromosomeIx][side], positions[chromosomeIx][side], -1),
                getLengths(segStarts[chromosomeIx][side], segEnds[chromosomeIx][side]),
                getMeans(segStarts[chromosomeIx][side], segEnds[chromosomeIx][side],
                        lnSignalSums[chromosomeIx][side]),
                getVars(segStarts[chromosomeIx][side], segEnds[chromosomeIx][side],
                        lnSignalSums[chromosomeIx][side], lnSignalSquareSums[chromosomeIx][side]));
        return arm;
    }

    private float getSum(int chromosomeIx, int arm, int segmentIx, float[][][] sums) {
        int start = segStarts[chromosomeIx][arm].get(segmentIx);
        int end = segEnds[chromosomeIx][arm].get(segmentIx);
        return sums[chromosomeIx][arm][end] - sums[chromosomeIx][arm][start];
    }

    public int getNumOfSignals(int chromosomeIx, int arm, int segmentIx) {
        return segEnds[chromosomeIx][arm].get(segmentIx) - segStarts[chromosomeIx][arm].get(segmentIx);
    }

    public int getLength(int chromosomeIx, int arm, int segmentIx) {
        return positions[chromosomeIx][arm][segEnds[chromosomeIx][arm].get(segmentIx) - 1] - //positions[chromosomeIx][arm][3149]
                positions[chromosomeIx][arm][segStarts[chromosomeIx][arm].get(segmentIx)];
    }

    public double getSegmentMean(int chromosomeIx, int arm, int segmentIx) {
        return ((double) getSum(chromosomeIx, arm, segmentIx, lnSignalSums)) /
                getNumOfSignals(chromosomeIx, arm, segmentIx) + shiftFactor;
    }

    public double getSegmentVar(int chromosomeIx, int arm, int segmentIx) {
        double mean = getSegmentMean(chromosomeIx, arm, segmentIx);
        float sos = getSum(chromosomeIx, arm, segmentIx, lnSignalSquareSums);
        return sos / getNumOfSignals(chromosomeIx, arm, segmentIx) - mean * mean;
    }

    private double[] getMeans(TIntList segStarts, TIntList segEnds, float[] sums) {
        double[] means = new double[segStarts.size()];
        for (int i = 0; i < segStarts.size(); ++i) {
            means[i] = (sums[segEnds.get(i)] - sums[segStarts.get(i)]) / (segEnds.get(i) - segStarts.get(i)) + shiftFactor;
        }
        return means;
    }

    private double[] getCorrectedMeans(TIntList segStarts, TIntList segEnds, float[] sums) {
        double[] means = getMeans(segStarts, segEnds, sums);
        for (int i = 0; i < segStarts.size(); ++i) {
            means[i] = correct(means[i]);
        }
        return means;
    }

    private double localIntegerCountFactor() {
        double factor = 0, n = 0;

        for (int chromosome = 1; chromosome <= 24; ++chromosome) {
            for (int arm = P; arm <= Q; ++arm) {
                if (segStarts[chromosome][arm] == null) {
                    continue;
                }
                int size = segStarts[chromosome][arm].size();
                for (int segment = 0; segment < size; ++segment) {
                    double mean = getSegmentMean(chromosome, arm, segment);
                    double correctedMean = correct(mean);
                    double delta = correctedMean - mean;
                    int nos = getNumOfSignals(chromosome, arm, segment);
                    factor += nos * delta;
                    n += nos;
                }
            }
        }
        return factor / n;
    }

    private boolean refineSegments() {
        boolean refined = false;

        for (int chromosome = 1; chromosome <= 24; ++chromosome) {
            for (int arm = P; arm <= Q; ++arm) {
                TIntList starts = segStarts[chromosome][arm];
                if (starts == null) {
                    continue;
                }
                TIntList ends = segEnds[chromosome][arm];
                double mean = getSegmentMean(chromosome, arm, 0);
                int prevSegCount = getCount(mean);
                for (int segment = 1; segment < starts.size(); ++segment) {
                    mean = getSegmentMean(chromosome, arm, segment);
                    int currSegCount = getCount(mean);
                    if (currSegCount == prevSegCount) {
                        // removing the current segment
                        starts.removeAt(segment);
                        ends.removeAt(segment - 1);
                        refined = true;
                    } else {
                        //						for (int i = 1; i <= MAX_REFINEMENT_WINDOW && e; ++i){

                        //						}
                    }
                    prevSegCount = currSegCount;
                }
            }
        }
        return refined;
    }

    public void printCountHistogram() {
        TIntList histogram = new TIntArrayList();
        int n = 0;
        float[][][] signalSums = new float[25][2][];
        for (int chromosome = 1; chromosome <= 24; ++chromosome) {
            for (int arm = P; arm <= Q; ++arm) {
                TIntList starts = segStarts[chromosome][arm];
                if (starts == null) {
                    continue;
                }
                signalSums[chromosome][arm] = prefixSumOfExponents(lnSignals[chromosome][arm]);
                for (int segment = 0; segment < starts.size(); ++segment) {
                    int currSegCount = getCount(getSegmentMean(chromosome, arm, segment));
                    float expEverage = getSum(chromosome, arm, segment, signalSums) /
                            getNumOfSignals(chromosome, arm, segment);

                    while (histogram.size() <= currSegCount) {
                        histogram.add(0);
                    }
                    int length = getLength(chromosome, arm, segment);
                    histogram.set(currSegCount, histogram.get(currSegCount) + length);
                    n += length;
                }
            }
        }

        for (int i = 0; i < histogram.size(); ++i) {
            System.out.println(i + ": " + (((double) histogram.get(i)) / n));
        }
    }

    public void writeToFile(String path) throws IOException {
        ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(path + ".gz")));
        out.writeObject(this);
        out.close();
    }

}
