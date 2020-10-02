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
- Additional foldback read arguments are optional. They are "-f=[f1,f2,...,fk]" to set the foldback vector,
"-E=<error model class name>" to set the foldback count error model, and "-W=<min weight>" to set the minimum foldback
weight. They function similarly to count vectors, where foldback weight bounds are calculated separately from count
vector weight bounds
- If "-f" is included, you can also explore the neighborhood of combining adjacent groups of segments that are
within some distance of each other, specified by "-c=<min combine weight>", according to the count weight function

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

6. Getting all BFB strings corresponding to the given counts and foldback reads:
-s -a -f=[0,1] [1,2]

6. Getting all BFB strings with counts at least 0.9 with respect to the given counts, and foldback read counts at least
0.3 with respect to the given foldback read vector, according to the polynomial count and fbr model respectively.
Adjacent segments can be combined if they are within 0.001 of each other.
-s -a -c=0.001 -e=PolynomialCountErrorModel -E=PolynomialFbrErrorModel -w=0.9 -W=0.3 -f=[1,0,2] [3,3,4]

 */

import bfbf.distances.Distances;
import bfbf.palindromes.BFBPalindrome;
import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.ErrorModel;
import bfbf.weights.FbrWeights;
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

    private static final String FOLDBACK = "f";
    private static final String FBR_ERROR_MODEL = "E";
    private static final String FBR_MIN_WEIGHT = "W";

    private static final String COMBINE_WEIGHT = "c";

    private static final String DISTANCE = "d";

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

        Option combineWeightOption = new Option(COMBINE_WEIGHT, true, "Minimum weight to combine segments");
        weightOption.setValueSeparator('=');
        weightOption.setArgName("minCombineWeight");
        options.addOption(combineWeightOption);

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

        Option fbrVectorOption = new Option(FOLDBACK, true,
                "Foldback read vector");
        fbrVectorOption.setValueSeparator('=');
        fbrVectorOption.setArgName("foldbackVect");
        options.addOption(fbrVectorOption);

        Option fbrWeightOption = new Option(FBR_MIN_WEIGHT, true,
                "Minimum solution foldback read weight");
        fbrWeightOption.setValueSeparator('=');
        fbrWeightOption.setArgName("minFbrWeight");
        options.addOption(fbrWeightOption);

        Option fbrErrorOption = new Option(FBR_ERROR_MODEL, true, "Foldback reads error model");
        fbrErrorOption.setValueSeparator('=');
        fbrErrorOption.setArgName("fbrErrorModel");
        options.addOption(fbrErrorOption);

        Option distanceOption = new Option(DISTANCE, true, "Distance model to combine fbr + count weights");
        distanceOption.setValueSeparator('=');
        distanceOption.setArgName("weightDistanceModel");
        options.addOption(distanceOption);

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

            double minCombineWeight = Double.parseDouble(getOptValue(cmd, COMBINE_WEIGHT, "1"));

            String errorModelClassName = getOptValue(cmd, ERROR_MODEL, "NoErrorModel");

            Class<?> clazz = Class.forName(ErrorModel.class.getPackageName() + "." + errorModelClassName);
            Constructor<?> ctor = clazz.getConstructor();
            errorModel = (ErrorModel) ctor.newInstance();

            FbrWeights fw = null;
            ErrorModel fbrErrorModel = null;
            double minFbrWeight = Double.parseDouble(getOptValue(cmd, FBR_MIN_WEIGHT, "1"));
            String fbrErrorModelClassName = getOptValue(cmd, FBR_ERROR_MODEL, "NoErrorModel");
            String fbrVect = null;
            File fbrInputFile = null;
            Class<?> fbrClazz = Class.forName(ErrorModel.class.getPackageName() + "." + fbrErrorModelClassName);
            Constructor<?> fbrCtor = fbrClazz.getConstructor();
            fbrErrorModel = (ErrorModel) fbrCtor.newInstance();
            Distances weightDistanceModel = null;
            String weightDistanceModelClassName = getOptValue(cmd, DISTANCE, "Distances");

            Class<?> distanceClazz = Class.forName(Distances.class.getPackageName() + "." + weightDistanceModelClassName);
            Constructor<?> distanceCtor = distanceClazz.getConstructor();
            weightDistanceModel = (Distances) distanceCtor.newInstance();


            if (cmd.hasOption(FOLDBACK)) {
                try {
                    fbrVect = getOptValue(cmd, FOLDBACK, inputStr);
                    fbrInputFile = new File(fbrVect);
                    if (fbrInputFile.isFile()) {
                        fbrVect = fileContent(fbrInputFile).trim();
                    }
                    fw = new FbrWeights(fbrVect, fbrErrorModel, minFbrWeight);
                } catch (IllegalArgumentException e) {
                    if (fbrInputFile.isFile()) {
                        fbrVect = " in " + fbrInputFile.getAbsolutePath();
                    } else {
                        fbrVect = ": " + fbrVect;
                    }
                    System.out.println("Invalid input format" + fbrVect);
                    System.out.println("Run with -h argument for help.");
                    System.exit(0);
                }
            }
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

                if (cmd.hasOption(FOLDBACK)) {
                    out.print("Foldback read input: ");
                    if (fbrInputFile.isFile()) {
                        out.println(fbrInputFile.getAbsolutePath());
                    } else {
                        out.println(fbrVect);
                    }

                    if (fbrErrorModel != null) {
                        out.print("Foldback read error model: " + fbrErrorModel.toString());
                        if (!(fbrErrorModel instanceof NoErrorModel)) {
                            out.println(", minimum solution weight: " + minFbrWeight);
                        } else {
                            out.println();
                        }
                    }
                    out.println("Weight distance error model: " + weightDistanceModel.toString());
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
                        if (cmd.hasOption(FOLDBACK)) {
                            FbrSolution solution = new FbrSolution();
                            int[] tmpFbrCounts = fw.getCounts();
                            int[] tmpCounts = w.getCounts();
                            ArrayList<int[]>[] allLists = getAllDiffLengthLists(tmpCounts, tmpFbrCounts, minCombineWeight, errorModel);
                            for (int j = 0; j < allLists[0].size(); j++) {
                                tmpCounts = allLists[0].get(j);
                                tmpFbrCounts = allLists[1].get(j);
                                int[] fbrCounts = new int[tmpFbrCounts.length + 1];
                                int[] counts = new int[tmpCounts.length + 1];
                                fbrCounts[0] = 0;
                                counts[0] = 1;
                                for (int i = 0; i < tmpFbrCounts.length; i++) {
                                    fbrCounts[i + 1] = 2 * tmpFbrCounts[i];
                                    counts[i + 1] = tmpCounts[i];
                                }
                                fw.updateCounts(fbrCounts, 0, fbrErrorModel, minFbrWeight);
                                w.updateCounts(counts, errorModel, minWeight);
                                //
                                // System.out.println("input: " + Arrays.toString(tmpCounts) + ", " + Arrays.toString(tmpFbrCounts));
                            /*for (int j = 0; j < counts.length; j++){
                                System.out.println("counts[" + j + "] = " + counts[j]);
                                System.out.println("weight of count = " + fw.getWeight(j, counts[j]));
                            }*/
                                solution = Signature.heaviestBFBVector(w, fw, minLength, minWeight,
                                        minFbrWeight, weightDistanceModel);
                                if (solution != null) {
                                    out.print("Heaviest BFB vector's count weight: " + solution.getWeight() +
                                            ", fold-back weight: " + solution.getFbrWeight() +
                                            " (combined: " + weightDistanceModel.combineWeights(solution.getWeight(), solution.getFbrWeight()) + ")");
                                    out.println(", counts:");
                                    out.println(Arrays.toString(solution.counts));
                                    out.println("fold-backs:");
                                    out.println(Arrays.toString(solution.displayFbrCounts()));
                                    out.println();
                                    break;
                                }
                            /*int[] test = {0,0,0,9,3,3,7};
                            int[] output ={0,2,0,9,5,4,9};
                            FbrWeights testWeights = new FbrWeights(test);
                            double currWt = 1;
                            for (int i = 0; i < test.length; i++){
                                System.out.println(i + ":" + output[i] + " vs " + test[i]);
                                currWt *= testWeights.getWeight(i, output[i]);
                                System.out.println(i + ": " + currWt);
                                System.out.println(i + ": " + fbrErrorModel.weight(output[i], test[i]));
                            }
                            System.out.println(currWt);
                            System.out.println(Arrays.toString(fw.minCounts));*/
                            }
                        }
                        else {
                            Solution solution = Signature.heaviestBFBVector(w, minLength, minWeight);
                            if (solution != null) {
                                out.println(solution.toString());
                                out.println("Heaviest BFB vector's weight: " + solution.getWeight() + ", counts:");
                                out.println(Arrays.toString(solution.counts));
                            }
                        }
                    }
                    break;
                case STRING:
                    if (cmd.hasOption(FOLDBACK)){
                        int currNumStrings = 0;
                        int[] tmpFbrCounts = fw.getCounts();
                        int[] tmpCounts = w.getCounts();
                        ArrayList<int[]>[] allLists = getAllDiffLengthLists(tmpCounts, tmpFbrCounts, minCombineWeight, errorModel);
                        for (int j = 0; j < allLists[0].size(); j++) {
                            if (j > 0 && Arrays.equals(allLists[0].get(j-1), allLists[0].get(j)) &&
                                Arrays.equals(allLists[1].get(j-1), allLists[1].get(j))){
                                continue;
                            }
                            tmpCounts = allLists[0].get(j);
                            tmpFbrCounts = allLists[1].get(j);
                            int[] fbrCounts = new int[tmpFbrCounts.length];
                            int[] counts = new int[tmpCounts.length];
                            for (int i = 0; i < tmpFbrCounts.length; i++) {
                                fbrCounts[i] = 2 * tmpFbrCounts[i];
                                counts[i] = tmpCounts[i];
                            }
                            fw.updateCounts(fbrCounts, 0, fbrErrorModel, minFbrWeight);
                            w.updateCounts(counts, errorModel, minWeight);
                            if (cmd.hasOption(ALL)) {
                                if (j == allLists[0].size() - 1){
                                    fw.setNumLoopsIndex(-1);
                                }
                                currNumStrings = allBFBStrings(w, fw, minWeight, minFbrWeight, minLength, -1, currNumStrings);
                            } else {
                                fw.setNumLoopsIndex(-1);
                                currNumStrings = allBFBStrings(w, fw, minWeight, minFbrWeight, minLength, 1, currNumStrings);
                                if (currNumStrings > 0){
                                    break;
                                }
                            }
                            /*counts[i] += 1;
                            if (i > 0){
                                counts[i-1] -= 1;
                            }*/

                            /*for (int j = 0; j < counts.length; j++){
                                System.out.println("counts[" + j + "] = " + counts[j]);
                                System.out.println("weight of count = " + fw.getWeight(j, counts[j]));
                            }*/
//			}
                        }
                    }else {
                        if (cmd.hasOption(ALL)) {
                            allBFBStrings(w, minWeight, minLength);
                        } else {
                            allBFBStrings(w, minWeight, minLength, 1);
                        }
                    }
                    break;
                default: // DECISION
                    if (cmd.hasOption(FOLDBACK)) {
                        FbrSolution solution = new FbrSolution();
                        int[] tmpFbrCounts = fw.getCounts();
                        int[] tmpCounts = w.getCounts();
                        ArrayList<int[]>[] allLists = getAllDiffLengthLists(tmpCounts, tmpFbrCounts, minCombineWeight, errorModel);
                        for (int j = 0; j < allLists[0].size(); j++) {
                            tmpCounts = allLists[0].get(j);
                            tmpFbrCounts = allLists[1].get(j);
                            int[] fbrCounts = new int[tmpFbrCounts.length + 1];
                            int[] counts = new int[tmpCounts.length + 1];
                            fbrCounts[0] = 0;
                            counts[0] = 1;
                            for (int i = 0; i < tmpFbrCounts.length; i++) {
                                fbrCounts[i + 1] = 2 * tmpFbrCounts[i];
                                counts[i + 1] = tmpCounts[i];
                            }
                            fw.updateCounts(fbrCounts, 0, fbrErrorModel, minFbrWeight);
                            w.updateCounts(counts, errorModel, minWeight);
                            //
                            // System.out.println("input: " + Arrays.toString(tmpCounts) + ", " + Arrays.toString(tmpFbrCounts));
                            /*for (int j = 0; j < counts.length; j++){
                                System.out.println("counts[" + j + "] = " + counts[j]);
                                System.out.println("weight of count = " + fw.getWeight(j, counts[j]));
                            }*/
                            solution = Signature.heaviestBFBVector(w, fw, minLength, minWeight,
                                    minFbrWeight, weightDistanceModel);
                            if (solution != null) {
                                break;
                            }
                        }
                        out.println("Result: " + (solution != null));
                    }else {
                        out.println("Result: " + (Signature.heaviestBFBVector(w, minLength, minWeight) != null));
                    }
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

    public static boolean canCombineCounts(int n1, int n2, double avg, double minWeight, ErrorModel errorModel){
        return ((errorModel.weight(n1, (int) (avg + 0.5)) >= minWeight ||
                errorModel.weight(n1, (int) (avg - 0.5)) >= minWeight) &&
                (errorModel.weight((int) (avg + 0.5), n1) >= minWeight ||
                        errorModel.weight((int) (avg - 0.5), n1) >= minWeight) &&
                (errorModel.weight(n1, n2) >= minWeight &&
                errorModel.weight(n2, n1) >= minWeight));
    }

    public static ArrayList<int[]>[] getAllDiffLengthLists(int[] n, int[] f, double minDiffWeight, ErrorModel errorModel){
        ArrayList<double[]> currLists = new ArrayList<>();
        ArrayList<Integer> currGroups = new ArrayList<>();
        ArrayList<double[]> prevLists = new ArrayList<>();
        ArrayList<Integer> prevGroups = new ArrayList<>();
        ArrayList<int[]> currFbrs = new ArrayList<>();
        ArrayList<int[]> prevFbrs = new ArrayList<>();
        double[] first = new double[1];
        first[0] = n[0];
        int[] fbrs = new int[1];
        fbrs[0] = f[0];
        prevLists.add(first);
        prevGroups.add(1);
        prevFbrs.add(fbrs);
        for (int l = 1; l < n.length; l++){
            currLists.clear();
            currGroups.clear();
            currFbrs.clear();
            ArrayList<Integer> indicesToAdd = new ArrayList<Integer>();
            for (int i = 0; i < prevLists.size(); i++){
                double[] prevCounts = prevLists.get(i);
                double[] counts = new double[prevCounts.length + 1];
                System.arraycopy(prevCounts, 0, counts, 0, prevCounts.length);
                counts[prevCounts.length] = n[l];

                int[] prevFbr = prevFbrs.get(i);
                int[] fbr = new int[prevFbr.length + 1];
                System.arraycopy(prevFbr, 0, fbr, 0, prevFbr.length);
                fbr[prevFbr.length] = f[l];

                currLists.add(counts);
                currGroups.add(1);
                currFbrs.add(fbr);

                int prevGroupSize = prevGroups.get(i);
                double tempAvg = (prevGroupSize * prevCounts[prevCounts.length - 1] + n[l])/(prevGroupSize + 1.0);
                if (canCombineCounts(n[l], n[l-1], tempAvg, minDiffWeight, errorModel)){
                    indicesToAdd.add(i);
                }
            }
            for (int i = 0; i < indicesToAdd.size(); i++){
                double[] prevCounts = prevLists.get(i);
                int prevGroupSize = prevGroups.get(i);
                double tempAvg = (prevGroupSize * prevCounts[prevCounts.length - 1] + n[l])/(prevGroupSize + 1.0);
                prevCounts[prevCounts.length - 1] = tempAvg;

                int[] prevFbr = prevFbrs.get(i);
                prevFbr[prevFbr.length - 1] += f[l];


                currLists.add(prevCounts);
                currGroups.add(prevGroupSize + 1);
                currFbrs.add(prevFbr);

            }
            prevLists = new ArrayList<>(currLists);
            prevGroups = new ArrayList<>(currGroups);
            prevFbrs = new ArrayList<>(currFbrs);
        }
        ArrayList<int[]> outputLists = new ArrayList<>();
        for (int i = 0; i < currLists.size(); i++){
            double[] currCounts = currLists.get(i);
            int[] counts = new int[currCounts.length];
            for (int j = 0; j < currCounts.length; j++){
                counts[j] = (int) (currCounts[j] + 0.5);
            }
            outputLists.add(counts);
        }

        ArrayList<int[]>[] output = new ArrayList[2];
        output[0] = outputLists;
        output[1] = currFbrs;

        return output;
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



    public static int allBFBStrings(Weights w, FbrWeights fw, double minWeight, double minFbrWeight, int minLength) {
        return allBFBStrings(w, fw, minWeight, minFbrWeight, minLength, System.out, -1, 0);
    }


    public static int allBFBStrings(Weights w, FbrWeights fw, double minWeight, double minFbrWeight,
                                     int minLength, PrintStream stream) {
        return allBFBStrings(w, fw, minWeight, minFbrWeight, minLength, stream, -1, 0);
    }

    public static int allBFBStrings(Weights w, FbrWeights fw, double minWeight, double minFbrWeight,
                                     int minLength, int maxStrings) {
        return allBFBStrings(w, fw, minWeight, minFbrWeight, minLength, System.out, maxStrings, 0);
    }

    public static int allBFBStrings(Weights w, FbrWeights fw, double minWeight, double minFbrWeight,
                                    int minLength, int maxStrings, int currNumStrings) {
        return allBFBStrings(w, fw, minWeight, minFbrWeight, minLength, System.out, maxStrings, currNumStrings);
    }

    public static void allBFBStrings(Weights w, double minWeight, int minLength, PrintStream stream, int maxStrings) {
        AllBFBStringPrinter handler = new AllBFBStringPrinter(w, null, 0, w.length(), minWeight, 0,
                stream, maxStrings);
        handler.handle(new PalindromeCollection(), w.length() - 1, 1, 1, -1, -1);
        stream.println("Total number of strings: " + handler.numOfStrings());
    }

    public static int allBFBStrings(Weights w, FbrWeights fw, double minWeight, double minFbrWeight, int minLength,
                                     PrintStream stream, int maxStrings, int currNumStrings) {
        AllBFBStringPrinter handler = new AllBFBStringPrinter(w, fw, 0, w.length(), minWeight, minFbrWeight,
                stream, maxStrings);
        handler.handle(new PalindromeCollection(), w.length() - 1, 1, 1, -1, -1);
        currNumStrings += handler.numOfStrings();
        if (fw.getNumLoopsIndex() == -1 || fw.getNumLoopsIndex() == w.length() - 1) {
            stream.println("Total number of strings: " + currNumStrings);
        }
        return currNumStrings;
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
