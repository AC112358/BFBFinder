import bfbf.*;
import bfbf.palindromes.BFBPalindrome;
import bfbf.palindromes.Palindrome;
import bfbf.palindromes.PalindromeCollection;
import bfbf.weights.NoErrorModel;
import bfbf.weights.Weights;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Tmp {

    static String[] strings = {
            "CBBBBCDEEDCBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABCDE",
            "CBAABBAABCDEEDCBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABCDE",
            "CBAAAABBAAAABCDEEDCBAAAAAAAAAAAAAAAAAAAAAAAAAABCDE",
            "CBAAAAAABBAAAAAABCDEEDCBAAAAAAAAAAAAAAAAAAAAAABCDE",
            "CBAAAAAAAABBAAAAAAAABCDEEDCBAAAAAAAAAAAAAAAAAABCDE",
            "CBAABCDEEDCBAAAAAAAAAAAAAAAABBAAAAAAAAAAAAAAAABCDE",
            "CCDEEDCBAAAAAAAAAAAAAAAABBAABBAAAAAAAAAAAAAAAABCDE",
            "CBAAAAAABCDEEDCBAAAAAAAAAAAAAABBAAAAAAAAAAAAAABCDE",
            "CCDEEDCBAAAAAAAAAAAAAABBAAAAAABBAAAAAAAAAAAAAABCDE",
            "CBAAAAAAAAAABBAAAAAAAAAABCDEEDCBAAAAAAAAAAAAAABCDE",
            "CBAAAAAAAAAABCDEEDCBAAAAAAAAAAAABBAAAAAAAAAAAABCDE",
            "CCDEEDCBAAAAAAAAAAAABBAAAAAAAAAABBAAAAAAAAAAAABCDE",
            "CCDEEDCBAAAAAAAAAABBAAAAAAAAAAAAAABBAAAAAAAAAABCDE",
            "CBAAAAAAAAAAAAAABCDEEDCBAAAAAAAAAABBAAAAAAAAAABCDE",
            "CBAAAAAAAAAAAABBAAAAAAAAAAAABCDEEDCBAAAAAAAAAABCDE",
            "CCDEEDCBAAAAAAAABBAAAAAAAAAAAAAAAAAABBAAAAAAAABCDE",
            "CBAAAAAAAAAAAAAAAAAABCDEEDCBAAAAAAAABBAAAAAAAABCDE",
            "CCDEEDCBAAAAAABBAAAAAAAAAAAAAAAAAAAAAABBAAAAAABCDE",
            "CBAAAAAAAAAAAAAAAAAAAAAABCDEEDCBAAAAAABBAAAAAABCDE",
            "CBAAAAAAAAAAAAAABBAAAAAAAAAAAAAABCDEEDCBAAAAAABCDE",
            "CCDEEDCBAAAABBAAAAAAAAAAAAAAAAAAAAAAAAAABBAAAABCDE",
            "CBAAAAAAAAAAAAAAAAAAAAAAAAAABCDEEDCBAAAABBAAAABCDE",
            "CCDEEDCBAABBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBAABCDE",
            "CBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABCDEEDCBAABBAABCDE",
            "CBAAAAAAAAAAAAAAAABBAAAAAAAAAAAAAAAABCDEEDCBAABCDE"
    };

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        //		//		int[] NCI_H508_14 = {1, 5, 3, 15, 7, 13, 5, 11, 7, 13, 8, 14, 4, 14, 4, 16};
        //		int[] NCI_H508_14 = {1, 5, 3, 15, 5, 13, 6, 12, 8, 12, 8, 14, 4, 16};
        //		String[] st1 = allBFBStrings(NCI_H508_14);
        //		//		testProjection();


        testAllFoldings();
//
        printSigCardinalities(1000);

        allStrings();

        System.exit(0);

        int[][] allCounts = new int[][]{
                {5, 3, 14, 5, 12, 5, 11, 7, 12, 7, 14, 4, 14, 4, 14, 5, 3, 5, 3}, // NCI-H508 chr14
                {4, 3, 4, 6, 5, 6, 7, 6, 7, 6, 5, 7, 6, 7, 6, 7, 3, 6, 8, 6, 7, 5}, // SK-CO-1 chr7
                {7, 10, 5, 7, 6, 5, 7, 5, 8, 0, 8, 5}, // SK-CO-1 chr14
                {4, 3, 5, 3, 1, 10, 5, 3, 5, 6, 9, 8, 7, 3, 1, 2, 1, 3, 4, 5, 1}, // SK-CO-1 chr18
                {3, 5, 4, 5, 6, 5, 4, 5, 4, 8}, // UM-UC-3 chr20
        };

        String[] sampleNames = {"NCI-H508_chr14", "SK-CO-1_chr7", "SK-CO-1_chr14", "SK-CO-1_chr18", "UM-UC-3_chr20"};

        int[] ploidies = {3, 2, 3, 2, 4};

        //		int[] counts = {3, 3, 4, 6, 34}; // Marcus' example
        final int minLength = 7;
        double maxNormalizedError = 0.2;

        for (int sampleIx = 0; sampleIx < allCounts.length; ++sampleIx) {
            int[] counts = BFB_Algo.fixPloidy(allCounts[sampleIx], ploidies[sampleIx]);

            String outputPath = "data/CGP/" + sampleNames[sampleIx] + ".bfb";

            final BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));

            String indexStr = "";
            String segmentStr = "";
            String countStr = tabDelim(counts);
            String originCountStr = tabDelim(allCounts[sampleIx]);

            for (int i = 0; i < counts.length; ++i) {
                indexStr += i + "\t";
                segmentStr += (char) ('A' + i) + "\t";
            }

            bw.append("> Sample:\t" + sampleNames[sampleIx] + "\n");
            bw.append("> Min substring length: " + minLength +
                    ", error model: " + Env.errorModel.getClass().getCanonicalName() +
                    ", max error: " + maxNormalizedError + "\n\n");

            bw.append("Original counts:\t" + originCountStr).append("\n");
            bw.append("Ploidy corrected (" + ploidies[sampleIx] + "):\t" + countStr).append("\n");
            bw.append("Segment indices:\t" + indexStr).append("\n");
            bw.append("Segment annotations:\t" + segmentStr).append("\n");


            long timeStart = System.currentTimeMillis();
            System.out.println("Counts: " + Arrays.toString(counts) + ", length: " + counts.length + ", scope: A-" + (char) ('A' + counts.length - 1));
            List<int[]> allBFBSubstring = BFBCalculator.allMaximalBFBSubstring(counts, minLength, maxNormalizedError);
            System.out.println("Number of maximal BFB substrings: " + allBFBSubstring.size());

            bw.append("> Maximal BFB sub-vector neighbors (" + allBFBSubstring.size() + "):\n");
//			bw.append("  \t\t").append(countStr).append("\n");
            for (int i = 0; i < allBFBSubstring.size(); ++i) {
                bw.append("#" + i + "\t" + tabDelim(allBFBSubstring.get(i)) + "\n");
            }

            final Set<String>[][] projections = new Set[counts.length][];
            for (int i = 1; i < counts.length; ++i) {
                projections[i] = new Set[i];
                for (int j = 0; j < i; ++j) {
                    projections[i][j] = new HashSet<String>();
                }
            }

            final int[] toBox = {0};
            final int[] neighborsBox = {0};
            final int[] sizeBox = {0};
            final Set<String> palindromeStrs = new HashSet<String>();

            Function<BFBPalindrome> collectStringsAndProjections = new Function<BFBPalindrome>() {
                StringBuilder sb = new StringBuilder();

                @Override
                public void execute(BFBPalindrome palindrome) {
                    ++neighborsBox[0];
                    int[] seq = palindrome.seq();
                    int depth = palindrome.depth();

                    int halfLength = seq.length / 2;
                    palindromeStrs.add(palindrome.toString((char) ('A' + toBox[0] - depth)).substring(0, halfLength));

                    char firstChar, secondChar;

                    int firstLastValue = depth;
                    if (depth == minLength) {
                        firstLastValue = 2;
                    }


                    for (int first = depth; first >= firstLastValue; --first) {
                        firstChar = (char) ('A' + toBox[0] - first);
                        for (int second = first - 1; second > 0; --second) {
                            secondChar = (char) ('A' + toBox[0] - second);
                            sb.setLength(0);
                            for (int k = 0; k < halfLength; ++k) {
                                if (seq[k] == first) {
                                    sb.append(firstChar);
                                } else if (seq[k] == second) {
                                    sb.append(secondChar);
                                }
                            }
                            if (projections[toBox[0] - second][toBox[0] - first].add(sb.toString())) {
                                ++sizeBox[0];
                            }
                        }

                    }
                }
            };

            int[] trimed; // = new int[minLength+1];
            int from, to;
            int stringIx = 0;
            for (int i = 0; i < allBFBSubstring.size(); ++i) {
                int[] substring = allBFBSubstring.get(i);
                neighborsBox[0] = 0;
                sizeBox[0] = 0;
                palindromeStrs.clear();
                System.out.print(++stringIx);
                int addedStrings = 0;
                int neighborStrings = 0;
                long currTime = System.currentTimeMillis();
                for (from = 0; substring[from] == 0; ++from) ;
                for (to = substring.length; substring[to - 1] == 0; --to) ;
                System.out.print(": neighboring counts: " + Arrays.toString(Arrays.copyOfRange(substring, from, to)));
                System.out.print(", length: " + (to - from) + ", scope: " + (char) ('A' + from) + "-" + (char) ('A' + to - 1));
                //			for (int start = from; start<to-minLength; ++start){


//				for (; to - from >= minLength; --to){
//					toBox[0] = to;
//					Set<PalindromeCollection> solutions = PalindromeCollection.allFoldings(
//							foldings, foldingWeights, substring, from, to, minLength, collectStringsAndProjections);
//
//					//				for (int end = start+minLength; end <= to; ++end){
//					//					currLength += substring[end-1];
//					//					trimed = Env.borrowIntArray(end-start+1);
//					//					trimed[0] = 1;
//					//					System.arraycopy(substring, start, trimed, 1, end-start);
//					////					String[] allFoldings = allBFBStrings(trimed, (char) ('A' + start));
//					//
//					//					byte[][] allSequences = allBFBSequences(trimed, currLength, (byte) (start-1));
//					//					neighborStrings += allSequences.length;
//					//					//			for (String folding : allFoldings){
//					//					for (int i = end-1; i > start; --i){
//					//						for (int j = i-1; j >= from; --j){
//					//							int prevSize = projections[i][j].size();
//					//							int totalLength = substring[j] + substring[i];
//					//							projections[i][j].addAll(reduce(allSequences, totalLength, 'A', j, i).keySet());
//					//							addedStrings += projections[i][j].size() - prevSize;
//					//						}
//					//					}
//					//					//			}
//					//					Env.returnIntArray(trimed);
//					//				}
//					CompositPalindrome.clearCached();
//				}

                bw.append("\n> BFB strings for neighbor #" + i + " (" + palindromeStrs.size() + ") " + Arrays.toString(substring) + "\n");
                String[] strings = palindromeStrs.toArray(new String[palindromeStrs.size()]);
                Arrays.sort(strings);
                for (int j = 0; j < strings.length; ++j) {
                    bw.append(strings[j]).append("\n");
                }

                System.out.println(", BFB strings: " + neighborsBox[0] + ", added projections: " + sizeBox[0] + ", time (sec): " + (System.currentTimeMillis() - currTime) / 1000);
            }
            System.out.println();


            bw.append("\n> All pairwise projections:\n");

            int min = Integer.MAX_VALUE;
            for (int i = 1; i < counts.length; ++i) {
                for (int j = 0; j < i; ++j) {
                    if (!projections[i][j].isEmpty()) {
                        min = Math.min(min, projections[i][j].size());
                        double alternations = 0;
                        bw.append("\n> " + (char) ('A' + j) + "-" + (char) ('A' + i) + " (" + projections[i][j].size() + ")\n");
                        String[] strings = projections[i][j].toArray(new String[projections[i][j].size()]);
                        Arrays.sort(strings);
                        for (int r = 0; r < strings.length; ++r) {
                            bw.append(strings[r]).append("\n");
                            char currChar = strings[r].charAt(0);
                            for (int q = 1; q < strings[r].length(); ++q) {
                                if (strings[r].charAt(q) != currChar) {
                                    ++alternations;
                                    currChar = strings[r].charAt(q);
                                }
                            }
                        }
                        bw.append("> Average alternations: " + (alternations / strings.length) + "\n");
                    }
                }
            }

            bw.close();

            System.out.println("Min projection pair size: " + min);
            for (int i = 1; i < counts.length; ++i) {
                for (int j = 0; j < i; ++j) {
                    if (projections[i][j].size() == min) {
                        System.out.println("All projections over " + (char) ('A' + j) + " and " + (char) ('A' + i) + ":");
                        for (String s : projections[i][j]) {
                            System.out.println(s);
                        }
                        System.out.println();
                    }
                }
            }

            System.out.println("Total time (sec):" + (System.currentTimeMillis() - timeStart) / 1000);
        }
    }

    private static String tabDelim(int[] counts) {
        return Arrays.toString(counts).replaceAll("[\\[\\]\\s]", "").replaceAll(",", "\t");
    }

    private static void testAllFoldings() throws IOException {
        int[] counts = {3, 3, 4, 6, 34};
        //		int[] counts = {1, 9, 5, 3, 4};
//		int[] counts = {1, 3, 121, 127, 126, 112, 108, 106, 100, 118, 54, 52, 86, 74, 54, 52, 10, 20, 18, 20, 20, 2};
        Weights weights = new NoErrorModel().getWeights(counts, 1);

        final Charset charset = StandardCharsets.UTF_8;
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos, true, charset.name());
        BFB.allBFBStrings(weights, 1, counts.length, ps);
        String[] st1 = new String(baos.toByteArray(), charset).split("\n");
        st1 = Arrays.copyOf(st1, st1.length - 1);
        ps.close();
        baos.close();

        String[] st2 = reverse(strings, 'A', 'E');
        Arrays.sort(st1);
        Arrays.sort(st2);

        for (int j = 0; j < st1.length; ++j) {
            if (!st1[j].equals(st2[j])) {
                System.out.println(j + "\n" + st1[j] + "\n" + st2[j] + "\n");
            }
        }

        System.out.println("Done!");
    }

    public static String[] allBFBStrings(int[] counts, final char firstChar) {
        //		Set<PalindromeCollection> prevSolutions = new HashSet<PalindromeCollection>();
        //		prevSolutions.add(new PalindromeCollection());
        //		Set<PalindromeCollection> currSolutions = new HashSet<PalindromeCollection>();
        //		Set<PalindromeCollection> seen = new HashSet<PalindromeCollection>();
        ////		ConvexedPalindrome gamma = new ConvexedPalindrome();
        //
        //		for (int l = counts.length-1; l>=0; --l){
        //			currSolutions.clear();
        //			seen.clear();
        //			for (PalindromeCollection solution : prevSolutions){
        //				solution.allFoldings1(currSolutions, seen, counts[l]);
        //			}
        //			for (PalindromeCollection solution : currSolutions){
        //				solution.wrap();
        //			}
        //			Set<PalindromeCollection> tmp = prevSolutions;
        //			prevSolutions = currSolutions;
        //			currSolutions = tmp;
        //		}

        final Set<String> allStrings = new HashSet<String>();
        Function<Palindrome> f = new Function<Palindrome>() {

            @Override
            public void execute(Palindrome arg) {
                String string = arg.toString(firstChar);
                allStrings.add(string.substring(0, string.length() / 2));
            }
        };

        List<PalindromeCollection> foldings = new ArrayList<>();
        TDoubleList foldingWeights = new TDoubleArrayList();
        Weights weights = new NoErrorModel().getWeights(counts, 1);

//		Set<PalindromeCollection> solutions = PalindromeCollection.allFoldings(
//				foldings, foldingWeights, weights, 0, counts.length, counts.length, f);

        String[] st1 = allStrings.toArray(new String[allStrings.size()]);

        //		int i=0;
        //		for (PalindromeCollection solution : solutions){
        //			String string = solution.get(0).toString((char) (firstChar-1));
        //			st1[i] = string.substring(1, string.length()/2);
        //			++i;
        //		}
        //		String[] st1 = new String[solutions.size()];
        Arrays.sort(st1);
        return st1;
    }

    private static byte[][] allBFBSequences(int[] counts, int currLength, byte firstChar) {

        Set<PalindromeCollection> solutions = null; //PalindromeCollection.allFoldings(counts, null, null, null, null);

        byte[][] sequences = new byte[solutions.size()][currLength];
        int i = 0;
        int depth;
        for (PalindromeCollection solution : solutions) {
            Palindrome palindrome = solution.get(0);
            int[] seq = palindrome.seq();
            depth = palindrome.depth();
            //			sequences[i] = new byte[seq.length/2-1];
            for (int j = currLength; j > 0; --j) {
                sequences[i][j - 1] = (byte) (depth - seq[j] + firstChar);
            }
            ++i;
        }
        //		Arrays.sort(sequences);
        return sequences;
    }


    private static String[] reverse(String[] strings, char bottom, char top) {
        String[] res = new String[strings.length];
        for (int i = 0; i < strings.length; ++i) {
            String reversed = "";
            for (int j = strings[i].length() - 1; j >= 0; --j) {
                reversed += (char) (bottom + top - strings[i].charAt(j));
            }
            res[i] = reversed;
        }
        return res;
    }

    private static void testProjection() {

        Arrays.sort(strings);

        System.out.println("Original strings (" + strings.length + "):");
        for (int i = 0; i < strings.length; ++i) {
            System.out.println(i + ")\t" + strings[i]);
        }
        System.out.println();

        char lastChar = 0;
        for (int i = 0; i < strings[0].length(); ++i) {
            lastChar = (char) Math.max(lastChar, strings[0].charAt(i));
        }

        for (char first = 'A'; first < lastChar; ++first) {
            for (char second = (char) (first + 1); second <= lastChar; ++second) {
                Map<String, List<Integer>> reduced = reduce(strings, first, second);
                System.out.println("Reduced strings with respect to " + first
                        + " and " + second + " (" + reduced.size() + "/" +
                        strings.length + "):");
                for (Entry<String, List<Integer>> entry : reduced.entrySet()) {
                    System.out.println(entry.getKey() + ", " + entry.getValue());
                }
                System.out.println();
            }
        }

        Map<String, List<Integer>> reduced = reduce(strings, 'C', 'D', 'E');
        System.out.println("Reduced strings with respect to " + "C, D, and E (" + reduced.size() + "/" +
                strings.length + "):");
        for (Entry<String, List<Integer>> entry : reduced.entrySet()) {
            System.out.println(entry.getKey() + ", " + entry.getValue());
        }
        System.out.println();

        reduced = reduce(strings, 'B', 'C', 'D', 'E');
        System.out.println("Reduced strings with respect to " + "B, C, D, and E (" + reduced.size() + "/" +
                strings.length + "):");
        for (Entry<String, List<Integer>> entry : reduced.entrySet()) {
            System.out.println(entry.getKey() + ", " + entry.getValue());
        }

        //		int[] counts = {1, 9, 5, 3, 4};
        //		System.out.println(BFBCalculator.searchBFB(counts));
    }


    private static Map<String, List<Integer>> reduce(String[] strs, char... toKeep) {
        String toKeepPtrnStr = "[";
        for (char c : toKeep) {
            toKeepPtrnStr += c;
        }
        toKeepPtrnStr += "]+";

        Pattern ptrn = Pattern.compile(toKeepPtrnStr);
        Map<String, List<Integer>> reduced = new TreeMap<String, List<Integer>>();

        for (int i = 0; i < strs.length; ++i) {
            String curr = "";
            Matcher matcher = ptrn.matcher(strs[i]);
            while (matcher.find()) {
                curr += matcher.group();
            }
            List<Integer> ls = reduced.get(curr);
            if (ls == null) {
                ls = new ArrayList<Integer>();
                reduced.put(curr, ls);
            }
            ls.add(i);
        }

        return reduced;
    }

    private static Map<String, List<Integer>> reduce(byte[][] sequences, int totalLength, char shift, int... toKeep) {
        Map<String, List<Integer>> reduced = new TreeMap<String, List<Integer>>();

        char[] charArr = new char[totalLength];
        int lastPosition = sequences[0].length - 1;
        int lastToKeepPosition = toKeep.length - 1;

        for (int i = 0; i < sequences.length; ++i) {
            int position = totalLength - 1;
            for (int j = lastPosition; j >= 0; --j) {
                for (int k = lastToKeepPosition; k >= 0; --k) {
                    if (sequences[i][j] == toKeep[k]) {
                        charArr[position] = (char) (shift + toKeep[k]);
                        --position;
                        break;
                    }
                }
            }
            String curr = new String(charArr);
            List<Integer> ls = reduced.get(curr);
            if (ls == null) {
                ls = new ArrayList<Integer>();
                reduced.put(curr, ls);
            }
            ls.add(i);
        }


        return reduced;
    }

    public static void printSigCardinalities(int maxSize) {
        List<Integer> startsWithZero = new ArrayList<Integer>(maxSize + 1);
        List<Integer> startsWithPositive = new ArrayList<Integer>(maxSize + 1);
        List<Integer> firstNonzeroIsPositive = new ArrayList<Integer>(maxSize + 1);
        List<Integer> atMostOne = new ArrayList<Integer>(maxSize + 1);
        List<Integer> allSeries = new ArrayList<Integer>(maxSize + 1);

        startsWithZero.add(1);
        startsWithPositive.add(0);
        firstNonzeroIsPositive.add(0);
        atMostOne.add(1);
        allSeries.add(1);

        List<Long> all = Signature.numOfSignatureLs(maxSize);


        for (int i = 1; i <= maxSize; ++i) {
            allSeries.add(allSeries.get(i - 1) + allSeries.get(i / 2));
            if (i % 2 == 0) {
                startsWithZero.add(allSeries.get(i / 2));
                startsWithPositive.add((allSeries.get(i) - startsWithZero.get(i)) / 2);
                firstNonzeroIsPositive.add(startsWithPositive.get(i) + firstNonzeroIsPositive.get(i / 2));
                atMostOne.add(startsWithZero.get(i));
            } else {
                startsWithZero.add(0);
                startsWithPositive.add((allSeries.get(i)) / 2);
                firstNonzeroIsPositive.add(startsWithPositive.get(i));
                atMostOne.add(firstNonzeroIsPositive.get(i / 2));
            }
            System.out.println(i + ", " + allSeries.get(i) +
                    ", (" + Signature.numOfSignatures(i, all) + ", " + firstNonzeroIsPositive.get(i)
                    + "), (" + Signature.numOfSignaturesSmallerThanUnit(i, all) + ", " + atMostOne.get(i) + ")");
        }
    }

    public static void allStrings() {
        int[] counts = {5, 3, 15, 6, 12, 6, 10, 8, 12, 8};
        strings = allBFBStrings(counts, 'A');
        testProjection();
        System.out.println(strings.length);
    }

}
