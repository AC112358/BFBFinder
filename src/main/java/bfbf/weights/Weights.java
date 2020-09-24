package bfbf.weights;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Weights {


    private static final Pattern countVecPtrn = Pattern.compile(
            "\\[\\s*(\\d+(\\s*,\\s*\\d+)*)?\\s*\\]");
    ;

    private static final String doubleRegex = "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?";
    private static final String weightsLineRegex = "(\\d+)\\s*:\\s*((" + doubleRegex + ")(\\s*,\\s*" + doubleRegex + ")*)";
    private static final Pattern weightsLinePtrn = Pattern.compile(weightsLineRegex);

    private static final Pattern weightsPattern = Pattern.compile(
            "\\{\\s*(" + weightsLineRegex + "(\\s*[,\\n]\\s*" + weightsLineRegex + ")*" + ")?\\s*\\}"
    );
    private static final int WEIGHTS_GROUP = 2;
    private static final int COUNT_GROUP = 1;


    public int[] minCounts;
    public int[] heaviestCounts;
    public double[][] weights;
    //protected int[] minCounts;
    //protected int[] heaviestCounts;
    //protected double[][] weights;
    protected int[] counts;

    protected Weights(int length) {
        minCounts = new int[length];
        heaviestCounts = new int[length];
        weights = new double[length][1];
        counts = new int[length];
    }

    public Weights(int[] counts) {
        this(counts.length);
        for (int i = 0; i < counts.length; ++i) {
            minCounts[i] = counts[i];
            heaviestCounts[i] = counts[i];
            weights[i][0] = 1;
            this.counts[i] = counts[i];
        }
    }

    public void updateCounts(int[] counts, ErrorModel errorModel, double minWeight){
        Weights w = errorModel.getWeights(counts, minWeight);
        this.minCounts = w.minCounts;
        this.weights = w.weights;
        this.heaviestCounts = w.heaviestCounts;
    }

    public Weights(String inputStr) throws IllegalArgumentException {
        this(inputStr, new NoErrorModel(), 1);
    }

    public Weights(String inputStr, ErrorModel errorModel, double minWeight){
        this(inputStr, errorModel, minWeight, countVecPtrn);
    }

    public Weights(String inputStr, ErrorModel errorModel, double minWeight, Pattern countVecPtrn) throws IllegalArgumentException {
        inputStr = inputStr.trim();
        Matcher matcher;
        matcher = countVecPtrn.matcher(inputStr);
        matcher = countVecPtrn.matcher(inputStr);
        if (matcher.matches()) {
            // Input is of a count vector format
            String[] countStrs = inputStr.replaceAll("[\\[\\]]", "").split("[,\\s]\\s*");
            int[] counts = new int[countStrs.length];
            this.counts = new int[countStrs.length];
            for (int i = 0; i < countStrs.length; ++i) {
                try {
                    counts[i] = Integer.parseInt(countStrs[i].trim());
                    this.counts[i] = counts[i];
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Invalid counts format: " + inputStr);
                }
            }
            Weights w = errorModel.getWeights(counts, minWeight);
            this.minCounts = w.minCounts;
            this.weights = w.weights;
            this.heaviestCounts = w.heaviestCounts;
        } else {
            matcher = weightsPattern.matcher(inputStr);
            if (!matcher.matches()) {
                throw new IllegalArgumentException("Invalid weights format:\n" + inputStr);
            }

            TIntList minCounts = new TIntArrayList();
            ArrayList<double[]> weights = new ArrayList<>();

            String weightStr = null;
            matcher = weightsLinePtrn.matcher(inputStr);
            try {
                while (matcher.find()) {
                    String[] weightStrs = matcher.group(WEIGHTS_GROUP).split(",");
                    double[] currWeights = new double[weightStrs.length];
                    for (int i = 0; i < weightStrs.length; ++i) {
                        currWeights[i] = Double.parseDouble(weightStrs[i].trim());
                    }
                    minCounts.add(Integer.parseInt(matcher.group(COUNT_GROUP).trim()));
                    weights.add(currWeights);
                }
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException(
//						"The following line contains " +
//						"the illegal weight term " + weightStr + " (weights must " +
//						"be numbers greater than 0 and smaller than or " +
//						"equal to 1):\n"
//						+ line
                );
            }

            this.minCounts = minCounts.toArray();
            this.weights = weights.toArray(new double[minCounts.size()][]);
            heaviestCounts = new int[minCounts.size()];
            processWeights();
        }
    }



    public static void main(String[] args) {
        String str = "\n\n   3:0.5264,1.555  \n6:0.2,0.4,0.8";
        Weights w = new Weights(str);
    }

    public int[] getCounts(){
        return counts;
    }

    public int length() {
        return minCounts.length;
    }

    public int[] getHeaviestCounts() {
        return heaviestCounts;
    }

    public int getMinCount(int ix) {
        if (ix < 0) return 1;
        else return minCounts[ix];
    }

    public int getMaxCount(int ix) {
        if (ix < 0) return 1;
        else return minCounts[ix] + weights[ix].length - 1;
    }

    public double getWeight(int ix, int count) {
        if (ix < 0) {
            if (count == 1) return 1;
            else return 0;
        } else {
            int j = count - minCounts[ix];
            if (j < 0 || j >= weights[ix].length) {
                return 0;
            } else {
                return weights[ix][j];
            }
        }
    }

    public double getMaxWeight(int ix) {
        if (ix < 0) return 1;
        else return getWeight(ix, heaviestCounts[ix]);
    }

    public void setWeights(int ix, int minCount, double... countWeights) {
        minCounts[ix] = minCount;
        weights[ix] = countWeights;
        processWeights(ix);
    }

    private void processWeights(int ix) {
        int heavistIx = 0;
        for (int i = 1; i < weights[ix].length; ++i) {
            if (weights[ix][i] > weights[ix][heavistIx]) {
                heavistIx = i;
            }
        }

        double factor = weights[ix][heavistIx];
        for (int i = 0; i < weights[ix].length; ++i) {
            weights[ix][i] /= factor;
        }
        heaviestCounts[ix] = minCounts[ix] + heavistIx;
    }

    public void processWeights() {
        for (int ix = 0; ix < weights.length; ++ix) {
            processWeights(ix);
        }
    }

    public void setWeights(int[] minCounts, double[][] weights) {
        this.minCounts = minCounts;
        this.weights = weights;
        heaviestCounts = new int[minCounts.length];
        processWeights();
    }

    public double weight(int[] counts) {
        int first = 0, last = counts.length;
        for (; first < last && counts[first] == 0; ++first) ;
        for (; first < last && counts[last - 1] == 0; --last) ;

        double weight = 1;
        for (int i = first; i < last; ++i) {
            weight *= getWeight(i, counts[i]);
        }
        return weight;
    }

    public int getHeaviestCount(int l) {
        if (l < 0) return 1;
        else return heaviestCounts[l];
    }

    public int getMinCount(int l, double w) {
        if (l < 0) return 1;
        int c = heaviestCounts[l];
        for (; c > 1 && getWeight(l, c - 1) >= w; --c) ;
        return c;
    }

    public int getMaxCount(int l, double w) {
        if (l < 0) return 1;
        int c = heaviestCounts[l];
        for (; getWeight(l, c + 1) >= w; ++c) ;
        return c;
    }
}
