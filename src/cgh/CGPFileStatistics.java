package cgh;

import bfbf.Env;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class CGPFileStatistics {

    public static final int P = 0, Q = 1;
    public static final String ENTRY_PATTERN_STR = "([^\\s,]*)";
    public static final String SEPARATOR_PATTERN_STR = "\\s*,\\s*";
    public int[][][] segmentCoordinates;
    public int[][][] segmentCount;
    public int[][][] signalNum;
    public double[][][] signalSum;
    public double[][][] signalSquareSum;
    public String[][][] additional;
    public String[] additionalTitles;
    public String[][][] statistics;

    public CGPFileStatistics(String path, String chromosomeTitle, String segmentIdentifierTitle,
                             String probLocationTitle, String probIdTitle, String signalTitle,
                             BooleanFunction<String> isValidProbe, String... additionalTitles) throws IOException {
        this(getSegmentCoordinatess(path, chromosomeTitle, segmentIdentifierTitle,
                probLocationTitle, CGHFileHandler.centromerPositions), path,
                chromosomeTitle, probLocationTitle, probIdTitle, signalTitle,
                isValidProbe, additionalTitles);
    }

    private CGPFileStatistics(int[][][] segmentCoordinates, String path,
                              String chromosomeTitle, String probLocationTitle, String probIdTitle,
                              String signalTitle, BooleanFunction<String> isValidProbe, String... additionalTitles) throws IOException {

        if (additionalTitles == null) {
            additionalTitles = new String[0];
        }

        int chromosomeNum = segmentCoordinates.length - 1;

        this.segmentCoordinates = segmentCoordinates;
        this.additionalTitles = additionalTitles;
        signalNum = new int[chromosomeNum + 1][2][];
        signalSum = new double[chromosomeNum + 1][2][];
        signalSquareSum = new double[chromosomeNum + 1][2][];
        additional = new String[chromosomeNum + 1][][];
        for (int i = 1; i <= chromosomeNum; ++i) {
            for (int j = 0; j < 2; ++j) {
                signalNum[i][j] = new int[segmentCoordinates[i][j].length];
                signalSum[i][j] = new double[segmentCoordinates[i][j].length];
                signalSquareSum[i][j] = new double[segmentCoordinates[i][j].length];
                additional[i] = new String[segmentCoordinates[i].length][additionalTitles.length];
            }
        }

        BufferedReader br = Env.getBufferedReader(path);
        String line = br.readLine(); // titles line
        Map<String, Integer> titlesMap = Env.columnTitles2Ixs(line, SEPARATOR_PATTERN_STR);

        int chromosomeIx = titlesMap.get(chromosomeTitle);
        int probIdIx = titlesMap.get(probIdTitle);
        int locationIx = titlesMap.get(probLocationTitle);
        int signalIx = titlesMap.get(signalTitle);
        int[] additionalIxs = new int[additionalTitles.length];
        for (int i = 0; i < additionalTitles.length; ++i) {
            additionalIxs[i] = titlesMap.get(additionalTitles[i]);
        }

        Pattern p = Env.getPattern(Math.max(Env.max(chromosomeIx, probIdIx, locationIx, signalIx), Env.max(additionalIxs)), SEPARATOR_PATTERN_STR, ENTRY_PATTERN_STR);
        String currCromosomeStr = "0";
        int currChromosome = 0;
        boolean isNewSegment = false;

        int segmentIx = 0;
        line = br.readLine();
        while (line != null) {
            Matcher matcher = p.matcher(line);
            if (matcher.find() && isValidProbe.execute(matcher.group(probIdIx))) {
                if (!currCromosomeStr.equals(matcher.group(chromosomeIx))) {
                    currCromosomeStr = matcher.group(chromosomeIx);
                    ++currChromosome;
                    segmentIx = 0;
                    isNewSegment = true;
                }
                int location = Integer.parseInt(matcher.group(locationIx));
                while (location >= segmentCoordinates[currChromosome][segmentIx]) {
                    ++segmentIx;
                    isNewSegment = true;
                }
                double signal = Double.parseDouble(matcher.group(signalIx));
                ++signalNum[currChromosome][segmentIx];
                signalSum[currChromosome][segmentIx] += signal;
                signalSquareSum[currChromosome][segmentIx] += signal * signal;
                if (isNewSegment) {
                    for (int i = 0; i < additionalIxs.length; ++i) {
                        additional[currChromosome][segmentIx][i] = matcher.group(additionalIxs[i]);
                    }
                    isNewSegment = false;
                }
            }
            line = br.readLine();
        }

    }

    private CGPFileStatistics(int[][] segmentCoordinates, String path,
                              String chromosomeTitle, String probLocationTitle, String probIdTitle,
                              BooleanFunction<String> isValidProbe, Accumulator<String, ?>... accumulators) throws IOException {


        int chromosomeNum = segmentCoordinates.length - 1;

        this.segmentCoordinates = segmentCoordinates;
        statistics = new String[chromosomeNum + 1][accumulators.length][];

        for (int i = 1; i <= chromosomeNum; ++i) {
            for (int j = 0; j < accumulators.length; ++j) {
                statistics[i][j] = new String[segmentCoordinates[i].length];
            }
        }

        BufferedReader br = Env.getBufferedReader(path);
        String line = br.readLine(); // titles line
        Map<String, Integer> titlesMap = Env.columnTitles2Ixs(line, SEPARATOR_PATTERN_STR);

        int chromosomeIx = titlesMap.get(chromosomeTitle);
        int probIdIx = titlesMap.get(probIdTitle);
        int locationIx = titlesMap.get(probLocationTitle);
        int[] titleIxs = new int[accumulators.length];
        for (int i = 0; i < accumulators.length; ++i) {
            titleIxs[i] = titlesMap.get(accumulators[i].inValueTitle);
        }

        Pattern p = Env.getPattern(Math.max(Env.max(chromosomeIx, probIdIx, locationIx), Env.max(titleIxs)), SEPARATOR_PATTERN_STR, ENTRY_PATTERN_STR);
        String currCromosomeStr = "0";
        int currChromosome = 0;

        int segmentIx = 0;
        line = br.readLine();
        while (line != null) {
            Matcher matcher = p.matcher(line);
            if (matcher.find() && isValidProbe.execute(matcher.group(probIdIx))) {
                String chromosomeStr = matcher.group(chromosomeIx);
                int location = Integer.parseInt(matcher.group(locationIx));
                if (!currCromosomeStr.equals(chromosomeStr) ||
                        location >= segmentCoordinates[currChromosome][segmentIx]) {
                    //new segment
                    for (int i = 0; i < accumulators.length; ++i) {
                        statistics[currChromosome][i][segmentIx] = accumulators[i].accumulatedValue().toString();
                        accumulators[i].clear();
                    }

                    if (!currCromosomeStr.equals(chromosomeStr)) {
                        currCromosomeStr = chromosomeStr;
                        ++currChromosome;
                        segmentIx = 0;
                    }
                    while (location >= segmentCoordinates[currChromosome][segmentIx]) {
                        ++segmentIx;
                    }

                }

                for (int i = 0; i < accumulators.length; ++i) {
                    accumulators[i].accumulate(matcher.group(titleIxs[i]));
                }
            }
            line = br.readLine();
        }

    }

    public static CGPFileStatistics make(String path, final String validProbIdentifier) throws IOException {
        return make(path, new BooleanFunction<String>() {
            @Override
            public boolean execute(String arg) {
                return arg.contains(validProbIdentifier);
            }
        });
    }

    public static CGPFileStatistics make(String path) throws IOException {
        return make(path, new BooleanFunction<String>() {
            @Override
            public boolean execute(String arg) {
                return true;
            }
        });
    }

    public static CGPFileStatistics make(String path, BooleanFunction<String> isValidProbe) throws IOException {
        return new CGPFileStatistics(path, "#CHR", "SEGMENT_COPYNUMBER", "POSITION",
                "AFFY_ID", "COPYNUMBER_INTENSITY", isValidProbe, "SEGMENT_COPYNUMBER");
    }

    private static int[][][] getSegmentCoordinatess(String path, String chromosomeTitle,
                                                    String segmentIdentifierTitle, String probLocationTitle,
                                                    int[] centromerPositions) throws IOException {

        BufferedReader br = Env.getBufferedReader(path);
        String line = br.readLine(); // titles line
        Map<String, Integer> titlesMap = Env.columnTitles2Ixs(line, SEPARATOR_PATTERN_STR);

        int chromosomeIx = titlesMap.get(chromosomeTitle);
        int locationIx = titlesMap.get(probLocationTitle);
        int segmentIdnIx = titlesMap.get(segmentIdentifierTitle);

        Pattern p = Env.getPattern(Env.max(chromosomeIx, locationIx, segmentIdnIx), SEPARATOR_PATTERN_STR, ENTRY_PATTERN_STR);

        List<TIntList[]> coordinates = new ArrayList<>();

        int segmentIx = 0;
        line = br.readLine();
        Matcher matcher = p.matcher(line);
        while (!matcher.find()) {
            line = br.readLine();
            matcher = p.matcher(line);
        }

        String currChromosomeStr = matcher.group(chromosomeIx);
        int currChromosomeIx = Integer.parseInt(currChromosomeStr);
        while (currChromosomeIx > coordinates.size()) {
            coordinates.add(null);
        }
        TIntList[] chromosomeCoordinates = new TIntArrayList[]{new TIntArrayList(), new TIntArrayList()};
        coordinates.add(chromosomeCoordinates);

        String prevCount = matcher.group(segmentIdnIx), currCount;
        int prevPosition = Integer.parseInt(matcher.group(locationIx)), currPosition;
        int arm = P;
        if (prevPosition >= centromerPositions[currChromosomeIx]) {
            arm = Q;
        }

        boolean isNewArm = true;
        chromosomeCoordinates[arm].add(prevPosition);

        line = br.readLine();
        while (line != null) {
            matcher = p.matcher(line);
            if (matcher.find()) {
                currPosition = Integer.parseInt(matcher.group(locationIx));
                currCount = matcher.group(segmentIdnIx);

                if (!currChromosomeStr.equals(matcher.group(chromosomeIx))) {
                    // A new chromosome
                    // Finalizing the previous chromosome segment coordinates:
                    chromosomeCoordinates[arm].add(prevPosition + 1);

                    // Setting the new chromosome parameters:
                    currChromosomeStr = matcher.group(chromosomeIx);
                    currChromosomeIx = Integer.parseInt(currChromosomeStr);
                    while (currChromosomeIx > coordinates.size()) {
                        coordinates.add(null);
                    }
                    chromosomeCoordinates = new TIntArrayList[]{new TIntArrayList(), new TIntArrayList()};
                    coordinates.add(chromosomeCoordinates);
                    arm = P;
                    if (currPosition >= centromerPositions[currChromosomeIx]) {
                        arm = Q;
                    }
                    chromosomeCoordinates[arm].add(currPosition);
                } else if (arm == P && currPosition >= centromerPositions[currChromosomeIx]) {
                    // A new arm
                    chromosomeCoordinates[P].add(prevPosition + 1);
                    arm = Q;
                    chromosomeCoordinates[Q].add(currPosition);
                } else {
                    if (!prevCount.equals(currCount)) {
                        // A new segment
                        chromosomeCoordinates[arm].add(currPosition);
                    }
                }
                prevCount = currCount;
                prevPosition = currPosition;
            }
            line = br.readLine();
        }

        chromosomeCoordinates[arm].add(prevPosition + 1);

        int[][][] segmentCoordinates = new int[coordinates.size()][2][];
        for (int i = 0; i < coordinates.size(); ++i) {
            chromosomeCoordinates = coordinates.get(i);
            if (chromosomeCoordinates != null) {
                segmentCoordinates[i][P] = chromosomeCoordinates[P].toArray();
                segmentCoordinates[i][Q] = chromosomeCoordinates[Q].toArray();
            }
        }
        return segmentCoordinates;
    }

    public static void main(String[] args) throws IOException {
        String identifier = "genotypes.csv.gz"; //".genotypes.csv.gz";
        String inputDir = "data/CGP/";
        String outputDir = "data/CGP/CN/DnaCopy/";
        File dir = new File(inputDir);
        String[] files = dir.list();
        for (String file : files) {
            if (file.contains(identifier)) {
                String path = inputDir + file;
                System.out.print("Collecting statistics from " + path);
                long start = System.currentTimeMillis();
                CGPFileStatistics statistics = make(path, "CN");
                statistics.writeToFile(outputDir + file.substring(0, file.indexOf(".") + 1) + "sgm");
                System.out.println(" (" + (System.currentTimeMillis() - start) / 1000 + " sec)");
            }
        }
    }

    public void writeToFile(String path) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(path));
        bw.append("Chromosome\tStart\tEnd\tSignals\tSumOfSignals\tSumOfSquareSignals");
        for (int i = 0; i < additionalTitles.length; ++i) {
            bw.append("\t").append(additionalTitles[i]);
        }
        bw.append("\n");

        int chromosomes = signalNum.length - 1;
        for (int i = 1; i <= chromosomes; ++i) {
            int coordinate = 0;
            int segments = signalNum[i].length;
            for (int j = 0; j < segments; ++j) {
                bw.append("" + i);
                bw.append("\t").append("" + coordinate);
                coordinate = segmentCoordinates[i][j];
                bw.append("\t").append("" + coordinate);
                bw.append("\t").append("" + signalNum[i][j]);
                bw.append("\t").append("" + signalSum[i][j]);
                bw.append("\t").append("" + signalSquareSum[i][j]);
                for (int k = 0; k < additional[i][j].length; ++k) {
                    bw.append("\t").append(additional[i][j][k]);
                }
                bw.append("\n");
            }
        }
        bw.close();
    }

}
