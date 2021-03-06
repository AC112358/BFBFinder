package bfbf;

/*

Main commandline arguments:
- Last argument should be a count vector of the format "[c1,c2,...,ck]" (c_i's are integers, no spaces).
- Strings describing weights or file paths (containing either count or weight strings) will be supported soon...
- Modes of run are either "Decision", "Vector" (finding neighboring BFB vectors), or "String".
If "-v" or "-s" arguments are supplied, the mode is set to "Vector" or "String" respectively, otherwise
it is "Decision".
- If "-a" argument is supplied, all solutions are returned, and otherwise a single solution is returned
("Decision" mode ignores this argument).
- To set error bounds, use "-e=<error model class name>" and "-w=<min weight>". Error model class names are
NoErrorModel (default), PoissonErrorModel, MaxRelativeErrorModel, and CanberraErrorModel. The min weight value
should be set between 0 (don't do that!) and 1. In short, a value of 1 allows no errors, while the closer
the value is to 0 more errors are allowed. Good chances there are still some bugs in this mechanism.

Usage:

Best way would be to use eclipse or intellij, and add the following 3 jars as dependencies:
BFBFinder/external/commons-cli-1.4-bin/commons-cli-1.4/commons-cli-1.4.jar
BFBFinder/external/commons-math3-3.6.1-bin/commons-math3-3.6.1/commons-math3-3.6.1.jar
BFBFinder/external/trove-3.0.3/3.0.3/lib/trove-3.0.3.jar

If you run from command line, set the working directory to "BFBFinder" (the project's root) and
use the following command, replacing "<arg list>" with concrete arguments:

java -classpath out/production/BFBFinder:external/commons-cli-1.4-bin/commons-cli-1.4/commons-cli-1.4.jar:external/commons-math3-3.6.1-bin/commons-math3-3.6.1/commons-math3-3.6.1.jar:external/trove-3.0.3/3.0.3/lib/trove-3.0.3.jar bfbf.BFB <arg list>

Some examples of possible arguments:

1. Valid BFB counts:
[3,3,4,6,34]

2. Invalid BFB counts:
[2,3,4,6,34]

3. Getting a single BFB string corresponding to the given counts:
-s [3,3,4,6,34]

3. Getting all BFB strings corresponding to the given counts:
-s -a [3,3,4,6,34]

4. Getting all BFB counts with weight at least 0.858 with respect to the given
counts, according to the Poisson error model:
-v -a -e=PoissonErrorModel -w=0.858 [4,3,4,6,34]

5. Getting all BFB strings corresponding to counts with weight at least 0.858 with respect to the given
counts, according to the Poisson error model:
-s -a -e=PoissonErrorModel -w=0.858 [4,3,4,6,34]

 */

import bfbf.palindromes.BFBPalindrome;
import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.ErrorModel;
import bfbf.weights.NoErrorModel;
import bfbf.weights.Weights;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import java.io.*;
import java.lang.reflect.Constructor;
import java.util.*;
import java.util.regex.Pattern;

public class BFB {

    //	private static final String DECISION = "t";
    private static final String VECTOR = "v";
    private static final String STRING = "s";

    private static final String ERROR_MODEL = "e";
    private static final String MIN_WEIGHT = "w";
    private static final String ALL = "a";
    private static final String SUBSTRING = "l";

    //	private static final String INPUT = "in";
    //	private static final String OUTPUT = "out";

    private static final String HELP = "h";
    private static final String VERBOSE = "b";

    private static final Weights[] weightsBox = new Weights[1];
    private static final Comparator<? super int[]> weightComparator = new Comparator<>() {

        @Override
        public int compare(int[] first, int[] second) {
            double firstWeight = weightsBox[0].weight(first);
            double secondWeight = weightsBox[0].weight(second);
            if (firstWeight > secondWeight) {
                return -1;
            } else if (firstWeight < secondWeight) {
                return 1;
            } else return 0;
        }
    };
    private static Pattern countVecPtrn = Pattern.compile("\\[?(\\d+(,\\d+)*)\\]?");

    ;

    /**
     * @param args
     */
    public static void main(String[] args) {
        // parse & set command line arguments
        Options options = new Options();

        //		options.addOption(DECISION, false, "Solves the decision problem variant");
        options.addOption(VECTOR, false, "Solves the vector search problem variant");
        options.addOption(STRING, false, "Solves the string search problem variant");
        options.addOption(ALL, false, "Reports all solutions");
        options.addOption(VERBOSE, false, "Verbose mode");

        //		options.addOption(INPUT, true, "A path to an input file or a count vector");
        //		options.addOption(OUTPUT, true, "A path to an output file");
        //
        Option weightOption = new Option(MIN_WEIGHT, true, "Minimum solution weight");
        weightOption.setValueSeparator('=');
        weightOption.setArgName("minWeight");
        options.addOption(weightOption);

        Option errorOption = new Option(ERROR_MODEL, true, "Error model");
        errorOption.setValueSeparator('=');
        errorOption.setArgName("errorModel");
        options.addOption(errorOption);

        Option substringOption = new Option(SUBSTRING, true, "Minimum number of " +
                "consecutive segments that can be reported as a BFB region (default " +
                "is all segments)");
        substringOption.setValueSeparator('=');
        substringOption.setArgName("minLength");
        options.addOption(substringOption);

        options.addOption(HELP, false, "Show help");

        Variant variant = Variant.DECISION;
        int minLength = 0;
        PrintStream out = System.out;
        boolean validArgs = true;
        Weights w = null;

        try {
            CommandLine cmd = new PosixParser().parse(options, args);
            if (cmd.hasOption(HELP) || args.length == 0) {
                printHelpAndExit();
            }

            ErrorModel errorModel = null;
            String inputStr, outputStr;

            if (args.length > 1 && !args[args.length - 2].startsWith("-")) {
                inputStr = args[args.length - 2];
                outputStr = args[args.length - 1];
                File outputFile = new File(outputStr);
                outputFile.createNewFile();

                if (!outputFile.isFile()) {
                    throw new Exception(outputStr + " is not a valid output file name.");
                } else if (!outputFile.canWrite()) {
                    throw new Exception("Cannot write output to " + outputStr + ".");
                } else {
                    out = new PrintStream(outputFile);
                }
            } else {
                inputStr = args[args.length - 1];
            }

            File inputFile = new File(inputStr);
            if (inputFile.isFile()) {
                inputStr = fileContent(inputFile).trim();
            }

//			w = new Weights(inputStr);

            double minWeight = Double.parseDouble(getOptValue(cmd, MIN_WEIGHT, "1"));

            String errorModelClassName = getOptValue(cmd, ERROR_MODEL, "NoErrorModel");

            Class<?> clazz = Class.forName(ErrorModel.class.getPackageName() + "." + errorModelClassName);
            Constructor<?> ctor = clazz.getConstructor();
            errorModel = (ErrorModel) ctor.newInstance();

//			Matcher matcher = countVecPtrn.matcher(inputStr);
//			if (matcher.matches()){
//				String[] countsStr = matcher.group(1).split(",");
//				int[] counts = new int[countsStr.length];
//				for (int i=0; i<counts.length; ++i){
//					counts[i] = Integer.parseInt(countsStr[i]);
//				}
//				w = new Weights(counts);
//			}
//			else {
            try {
                w = new Weights(inputStr, errorModel, minWeight);
            } catch (IllegalArgumentException e) {
                if (inputFile.isFile()) {
                    inputStr = " in " + inputFile.getAbsolutePath();
                } else {
                    inputStr = ": " + inputStr;
                }
                System.out.println("Invalid input format" + inputStr);
                System.out.println("Run with -h argument for help.");
                System.exit(0);
            }
//			}

            if (cmd.hasOption(VECTOR)) {
                variant = Variant.VECTOR;
            }
            if (cmd.hasOption(STRING)) {
                if (variant == Variant.VECTOR) {
                    System.out.println("Only one arguments among '" + VECTOR +
                            "' and '" + STRING + "' is allowed. Existing program...");
                }
                variant = Variant.STRING;
            }

            if (cmd.hasOption(SUBSTRING)) {
                minLength = Math.min(Integer.parseInt(getOptValue(cmd, SUBSTRING, "")), w.length());
            } else {
                minLength = w.length();
            }

            if (cmd.hasOption(VERBOSE)) {
                String modeStr = "";
                if (cmd.hasOption(ALL)) {
                    modeStr = ", exhaustive mode";
                }

                String variantStr;
                if (variant == Variant.VECTOR) {
                    variantStr = "vector search" + modeStr;
                } else if (variant == Variant.STRING) {
                    variantStr = "string search" + modeStr;
                } else {
                    variantStr = "decision";
                }
                out.println("BFB variant: " + variantStr);

                out.print("Input: ");
                if (inputFile.isFile()) {
                    out.println(inputFile.getAbsolutePath());
                } else {
                    out.println(inputStr);
                }

                out.print("Number of input segments: " + w.length());

                if (minLength == w.length()) {
                    out.println(", a solution must include all input segments");
                } else {
                    out.println(", minimum number of consecutive segments in a solution: " + minLength);
                }

                if (errorModel != null) {
                    out.print("Error model: " + errorModel.toString());
                    if (!(errorModel instanceof NoErrorModel)) {
                        out.println(", minimum solution weight: " + minWeight);
                    } else {
                        out.println();
                    }
                }
                out.println();
            }

            switch (variant) {
                case VECTOR:
                    if (cmd.hasOption(ALL)) {
                        List<int[]> allHeavyBFBSubVectors = Signature.allHeavyBFBSubVectors(w, minWeight, minLength);
                        if (allHeavyBFBSubVectors.isEmpty()) {
                            out.println("No solution found");
                        } else {
                            weightsBox[0] = w;
                            Collections.sort(allHeavyBFBSubVectors, weightComparator);
                            out.println("Solutions (" + allHeavyBFBSubVectors.size() + "):");
                            for (int[] counts : allHeavyBFBSubVectors) {
                                out.println(Arrays.toString(counts));
                            }
                        }

                    } else {
                        Solution1 solution = Signature.heaviestBFBVector(w, minLength, minWeight);
                        if (solution != null) {
                            out.println("Heaviest BFB vector's weight: " + solution.getWeight() + ", counts:");
                            out.println(Arrays.toString(solution.counts));
                        }
                    }
                    break;
                case STRING:
                    if (cmd.hasOption(ALL)) {
                        allBFBStrings(w, minWeight, minLength);
                    } else {
                        allBFBStrings(w, minWeight, minLength, 1);
                    }
                default: // DECISION
                    out.println("Result: " + (Signature.heaviestBFBVector(w, minLength, minWeight) != null));
            }

        }
//		catch (ParseException e) {
//			validArgs = false;
//			System.out.println("Cannot parse arguments " + Arrays.toString(args));
//			e.printStackTrace();
//		} catch (IOException e) {
//			validArgs = false;
//			System.out.println(e.getMessage());
//			e.printStackTrace();
//		}
        catch (Exception e) {
            validArgs = false;
            System.out.println(e.getMessage());
//			e.printStackTrace();
        }

        if (!validArgs) {
            System.out.println("Run with -h argument for help.");
            System.exit(0);
        }

    }

    private static String getOptValue(CommandLine cmd, String optName, String defaultVal) {
        return cmd.getOptionValue(optName, defaultVal).replaceAll("=", "");
    }

    private static String fileContent(File inputFile) throws IOException {
        StringBuilder sb = new StringBuilder();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = br.readLine();
        while (line != null) {
            sb.append(line).append("\n");
            line = br.readLine();
        }
        return sb.toString();
    }

    private static void printHelpAndExit() {
        String usage = "TODO: usage description";
        System.out.println("\n" + usage);
        System.exit(0);
    }

    public static boolean isBFB(String str) {
        str = str.toUpperCase();
        int k = 1;
        int length = str.length();
        for (; k < length && str.charAt(k) == str.charAt(k - 1) + 1; ++k) ;
        return isBFB(str, k, length);

    }

    private static boolean isBFB(String str, int k, int length) {
        if (k == length) {
            return true;
        } else {
            int min = Math.max(k, length / 2);
            for (int i = length - 1; i >= min; --i) {
                if (isPalindromeCenter(str, i, length)) {
                    return isBFB(str, k, i);
                }
            }
        }
        return false;
    }

    private static boolean isPalindromeCenter(String str, int i, int length) {
        for (int j = 0; i + j < length; ++j) {
            if (str.charAt(i + j) != str.charAt(i - j - 1)) {
                return false;
            }
        }
        return true;
    }

    public static List<BFBPalindrome> allBFBPalindormes(Weights w, double minWeight, int minLength) {
        List<SigCurve> oracles = SigCurve.sigCurves(w, minWeight);
        List<PalindromeCollection> prevCollections = new ArrayList<>();
        List<PalindromeCollection> currCollections = new ArrayList<>();
        List<PalindromeCollection> tmpCollections;
        currCollections.add(new PalindromeCollection());

        TDoubleList prevWeights = new TDoubleArrayList();
        TDoubleList currWeights = new TDoubleArrayList();
        TDoubleList tmpWeights;
        currWeights.add(1);

        for (int l = w.length() - 1; l >= 0; --l) {
            tmpCollections = currCollections;
            currCollections = prevCollections;
            prevCollections = tmpCollections;

            tmpWeights = currWeights;
            currWeights = prevWeights;
            prevWeights = tmpWeights;

            currCollections.clear();
            currWeights.clear();


            for (int i = 0; i < prevCollections.size(); ++i) {
                PalindromeCollection collection = prevCollections.get(i);
                double weight = prevWeights.get(i);
                collection.allFoldings(currCollections, currWeights, w, l, weight, minWeight, oracles.get(l));
            }

            for (PalindromeCollection collection : currCollections) {
                collection.wrap();
            }
        }

        prevCollections.clear();

        for (int i = 0; i < currCollections.size(); ++i) {
            PalindromeCollection collection = currCollections.get(i);
            double weight = currWeights.get(i);
            collection.allFoldings(prevCollections, currWeights, w, -1, weight, minWeight, oracles.get(0));
        }


        List<BFBPalindrome> palindromes = new ArrayList<>();
        for (PalindromeCollection collection : prevCollections) {
            assert collection.size() == 1;
            palindromes.add(collection.get(0));
        }
        return palindromes;
    }

    public static void allBFBStrings(Weights w, double minWeight, int minLength) {
        allBFBStrings(w, minWeight, minLength, System.out, -1);
    }

    public static void allBFBStrings(Weights w, double minWeight, int minLength, PrintStream stream) {
        allBFBStrings(w, minWeight, minLength, stream, -1);
    }

    public static void allBFBStrings(Weights w, double minWeight, int minLength, int maxStrings) {
        allBFBStrings(w, minWeight, minLength, System.out, maxStrings);
    }

    public static void allBFBStrings(Weights w, double minWeight, int minLength, PrintStream stream, int maxStrings) {
        AllBFBStringPrinter handler = new AllBFBStringPrinter(w, 0, w.length(), minWeight, stream, maxStrings);
        handler.handle(new PalindromeCollection(), w.length() - 1, 1);
        stream.println("Total number of strings: " + handler.numOfStrings());
    }

    public static void allBFBPairwiseProjections(Weights w, double minWeight, int minLength) {
        allBFBPairwiseProjections(w, minWeight, minLength, System.out, -1);
    }

    public static void allBFBPairwiseProjections(Weights w, double minWeight, int minLength, int maxStrings) {
        allBFBPairwiseProjections(w, minWeight, minLength, System.out, maxStrings);
    }

    public static void allBFBPairwiseProjections(Weights w, double minWeight, int minLength, PrintStream stream) {
        allBFBPairwiseProjections(w, minWeight, minLength, stream, -1);

    }

    public static void allBFBPairwiseProjections(Weights w, double minWeight, int minLength, PrintStream stream, int maxStrings) {
        AllPairwiseBFBStringAccumulator handler = new AllPairwiseBFBStringAccumulator(w, 0, w.length(), minWeight, maxStrings);
        handler.handle(new PalindromeCollection(), w.length() - 1, 1);
        handler.printAll(stream);
    }

    public static int countAllBFBStrings(Weights w, double minWeight, int minLength) {
        return countAllBFBStrings(w, minWeight, minLength, -1);
    }

    public static int countAllBFBStrings(Weights w, double minWeight, int minLength, int maxCount) {
        int[] countBox = {0};
        FoldingHandler handler = new AllBFBStringCounter(w, 0, w.length(), minWeight, countBox, maxCount);
        handler.handle(new PalindromeCollection(), w.length() - 1, 1);
        return countBox[0];
    }

    private static enum Variant {DECISION, VECTOR, STRING}
}
