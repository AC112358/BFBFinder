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

import bfbf.weights.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * The main entrance point for running the package's algorithms from the
 * command line.
 * <p>
 * Usage: java bfb.BFB_Algo counts:[<a count vector, comma delimited>] mode:<activation mode>
 * model:<error model name> maxError:<max error value>
 * <p>
 * Example: java bfb.BFB_Algo counts:[10,5,3,4,7,4] mode:substring model:Poisson maxError:0.22
 * <p>
 * Remarks:
 * <br>
 * (1) All arguments are optional. If "counts" is unspecified, counts are taken from the prompt.
 * <br>
 * (2) Available modes:
 * <br> - "decision" (default) - prints "true" if and only if the input vector
 * is a BFB count vector;
 * <br> - "search" - prints a BFB string corresponding to the input vector if
 * it is a BFB count vector, and otherwise prints "FAILED";
 * <br> - "suffix" - prints the longest BFB vector which approximates a suffix
 * of the input vector;
 * <br> - "substring" - prints the longest BFB vector which approximates a
 * consecutive sub-vector of the input vector;
 * <br> - "distance" - prints the minimum distances from BFB count vectors for
 * all substrings of the input vector.
 * <br><br>
 * (3) Available error models: "Poisson" (default), "noError", "Canberra", and
 * "maxRelative". Default  maxError is 0.
 * <br>
 * (4) Both "decision" and "search" modes ignore error model and max error
 * arguments if specified. For modes "suffix" and "substring", a BFB vector is
 * said to approximates the input vector if the distance of the input from the
 * BFB vector is at most the maximum error with respect to the error model. For
 * the "distance" mode, the distancesare computed with respect to the specified
 * error model, and sub-vectors whose minimum distance from a BFB vector is
 * more than the specified maximum error are considered to have an infinite
 * distance from a BFB vector.
 *
 * @author Shay Zakov
 */
public class BFB_Algo {


    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        if (args.length == 1 && args[0].equals("-man")) {
            System.out.println(format());
            return;
        }

        Env.errorModel = getErrorModel(args);
        double maxError = getMaxError(args);
        List<int[]> countsList = getCounts(args);
        PrintStream outStream = getOut(args);
        int mode = getMode(args);
        int ploidy = getPloidy(args);

        int[] counts;

        if (countsList != null) {
            for (int i = 0; i < countsList.size(); ++i) {
                counts = countsList.get(i);
                runMode(maxError, outStream, mode, ploidy, counts, null);
            }
        } else {
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            counts = parseCounts(br.readLine());
            while (counts != null) {
                runMode(maxError, outStream, mode, ploidy, counts, null);
                counts = parseCounts(br.readLine());
            }
            System.out.println("Non-count string received. Quitting.");
        }
    }

    private static void runMode(double maxError, PrintStream outStream,
                                int mode, int ploidy, int[] counts, double[][] countWeights) {
        outStream.println("Input vector:            " + Arrays.toString(counts));
        if (ploidy != 1) {
            counts = fixPloidy(counts, ploidy);
            outStream.println("After ploidy correction: " + Arrays.toString(counts));
        }
        if (mode == Env.DECISION) {
            outStream.println(BFBCalculator.decisionBFB(counts));
        } else if (mode == Env.SEARCH) {
            outStream.println(BFBCalculator.searchBFB(counts));
        } else if (mode == Env.DISTANCE) {
            double[][] distances = BFBCalculator.distanceBFB(counts, countWeights, maxError);
            for (int j = 0; j < counts.length; ++j) {
                outStream.println(Arrays.toString(distances[j]));
            }
        } else {
            outStream.println(BFBCalculator.longestBFB(counts, countWeights, maxError, mode));
        }
        outStream.println();
    }

    private static int getPloidy(String[] args) {
        int i = argIx(args, "ploidy:");
        int ploidy = 1;

        if (i < args.length) {
            String ploidyStr = args[i].substring(7);
            try {
                ploidy = Integer.parseInt(ploidyStr);
            } catch (NumberFormatException e) {
                ploidy = -1;
            }
            if (ploidy < 1) {
                System.err.println("Illegal ploidy value: " + ploidyStr
                        + ". Setting ploidy to 1.");
                ploidy = 1;
            }
        }
        return ploidy;
    }

    private static int getMode(String[] args) {
        int i = argIx(args, "mode:");
        int mode = Env.SUFFIX;

        if (i < args.length) {
            String modeStr = args[i].substring(5);
            if (modeStr.equalsIgnoreCase("decision")) {
                mode = Env.DECISION;
            } else if (modeStr.equalsIgnoreCase("search")) {
                mode = Env.SEARCH;
            } else if (modeStr.equalsIgnoreCase("suffix")) {
                mode = Env.SUFFIX;
            } else if (modeStr.equalsIgnoreCase("substring")) {
                mode = Env.SUBSTRING;
            } else if (modeStr.equalsIgnoreCase("distance")) {
                mode = Env.DISTANCE;
            } else {
                System.err.println("Unknown mode: " + modeStr
                        + ". Setting mode to \"suffix\".");
                mode = Env.SUFFIX;
            }
        }
        return mode;
    }

    private static double getMaxError(String[] args) {
        int i = argIx(args, "maxError:");
        double maxError = 0;
        if (i < args.length) {
            try {
                maxError = Double.parseDouble(args[i].substring(9));
            } catch (NumberFormatException e) {
                System.err.println("Illegal maxError argument: " + args[i].substring(9)
                        + ". Setting maxError to 0.");
            }
            if (maxError < 0) {
                System.err.println("maxError has to be nonnegative."
                        + " Setting maxError to 0.");
                maxError = 0;
            }
        }
        return maxError;
    }

    private static ErrorModel getErrorModel(String[] args) {
        int i = argIx(args, "model:");
        ErrorModel errorModel;
        if (i == args.length) {
            errorModel = new PoissonErrorModel();
        } else {
            String modelName = args[i].substring(6);
            if (modelName.equalsIgnoreCase("Poisson")) {
                errorModel = new PoissonErrorModel();
            } else if (modelName.equalsIgnoreCase("noError")) {
                errorModel = new NoErrorModel();
            } else if (modelName.equalsIgnoreCase("maxRelative")) {
                errorModel = new MaxRelativeErrorModel();
            } else if (modelName.equalsIgnoreCase("Canberra")) {
                errorModel = new CanberraErrorModel();
            } else {
                System.err.println("Unknown error model name: " + modelName +
                        ". Using Poisson model.");
                errorModel = new PoissonErrorModel();
            }
        }
        return errorModel;
    }


    private static String format() {
        return "Usage: java bfb.BFB_Algo " +
                "counts:[<a count vector, comma delimited>] " +
                "mode:<activation mode> " +
                "model:<error model name> " +
                "maxError:<max error value>\n" +
                //
                "\nExample: java bfb.BFB_Algo " +
                "counts:[10,5,3,4,7,4] " +
                "mode:substring " +
                "model:Poisson " +
                "maxError:0.22 \n" +
                //
                "\nRemarks: " +
                "\n(1) All arguments are optional. If \"counts\" is " +
                "unspecified, counts are taken from the prompt. " +
                //
                "\n(2) Available modes: " +
                //
                "\n\t\"decision\" (default) - prints \"true\" if and " +
                "only if the input vector is a BFB count vector; " +
                //
                "\n\t\"search\" - prints a BFB string corresponding to the " +
                "input vector if it is a BFB count vector, and otherwise " +
                "prints \"FAILED\"; " +
                //
                "\n\t\"suffix\" - prints the longest BFB vector which " +
                "approximates a suffix of the input vector; " +
                //
                "\n\t\"substring\" - prints the longest BFB " +
                "vector which approximates a consecutive sub-vector " +
                "of the input vector; " +
                //
                "\n\t\"distance\" - prints the minimum distances from BFB " +
                "count vectors for all substrings of the input vector. " +
                //
                "\n\n(3) Available error models: \"Poisson\" (default), \"noError\", " +
                "\"Canberra\", and \"maxRelative\". Default  maxError is 0." +
                //
                "\n(4) Both \"decision\" and \"search\" modes ignore error " +
                "model and max error arguments if specified. For modes \"suffix\" " +
                "and \"substring\", a BFB vector is said " +
                "to approximates the input vector if the distance of the input " +
                "from the BFB vector is at most the maximum error with respect " +
                "to the error model. For the \"distance\" mode, the distances" +
                "are computed with respect to the specified error model, and " +
                "sub-vectors whose minimum distance from a BFB vector is more " +
                "than the specified maximum error are considered to have an " +
                "infinite distance from a BFB vector."
                ;
    }

    private static List<int[]> getCounts(String[] args) {
        int i = argIx(args, "counts");
        if (i < args.length) {
            List<int[]> countsList = new ArrayList<>();
            String countStr = args[i].substring(7, args[i].length());
            File file = new File(countStr);
            if (file.isFile()) {
                try {
                    BufferedReader br = new BufferedReader(new FileReader(file));
                    for (String line = br.readLine(); line != null; line = br.readLine()) {
                        int[] counts = parseCounts(line);
                        if (counts != null) {
                            countsList.add(counts);
                        } else {
                            System.err.println("Illegal count vector line: " + line + ". Skeeping line...");
                        }
                    }
                    br.close();
                } catch (Exception e) {
                    System.err.println("Error when trying to read "
                            + file.toString() + ":\n" + e.getMessage());
                    System.exit(1);
                }
            } else {
                countsList.add(parseCounts(countStr));
            }
            return countsList;
        } else return null;
    }

    private static int[] parseCounts(String countStr) {
        int i;
        String[] countStrs = countStr.replaceAll("[\\[\\]]", "").split("[,\\s]\\s*");
        int[] counts = new int[countStrs.length];
        for (i = 0; i < countStrs.length; ++i) {
            try {
                counts[i] = Integer.parseInt(countStrs[i]);
            } catch (NumberFormatException e) {
                return null;
            }
        }
        return counts;
    }

    @SuppressWarnings("resource")
    private static PrintStream getOut(String[] args) {
        PrintStream outStream = System.out;
        int i = argIx(args, "out");
        if (i < args.length) {
            try {
                outStream = new PrintStream(args[i].substring(4));
            } catch (IOException e) {
                e.printStackTrace();
                outStream = System.out;
            }
        }
        return outStream;
    }

    private static int argIx(String[] args, String prefix) {
        int i = 0;
        while (i < args.length && !args[i].startsWith(prefix)) {
            ++i;
        }
        return i;
    }

    public static int[] fixPloidy(int[] counts, int ploidy) {
        if (ploidy > 1) {
            counts = Arrays.copyOf(counts, counts.length);
            for (int i = 0; i < counts.length; ++i) {
                counts[i] -= (ploidy - 1);
                if (counts[i] < 0) {
                    counts[i] = 0;
                }
            }
        }
        return counts;
    }

    public static void divideByMaxWeight(double[][] countWeights) {
        int maxIx;
        for (int i = 0; i < countWeights.length; ++i) {
            maxIx = 0;
            for (int j = 1; j < countWeights[i].length; ++j) {
                if (countWeights[i][maxIx] < countWeights[i][j]) {
                    maxIx = j;
                }
            }
            for (int j = 0; j < countWeights[i].length; ++j) {
                countWeights[i][j] /= countWeights[i][maxIx];
            }
        }
    }


}
