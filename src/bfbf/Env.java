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

//import cgh.CGHFileHandler;

import bfbf.weights.ErrorModel;
import bfbf.weights.PoissonErrorModel;
import bfbf.weights.Weights;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * An auxiliary class.
 *
 * @author Shay Zakov
 */
public class Env {

    public static final int SUFFIX = 0;
    public static final int SUBSTRING = 1;
    public static final int DISTANCE = 2;
    public static final int DECISION = 3;
    public static final int SEARCH = 4;

    public static final Map<String, String> chrStr2Int = new HashMap<String, String>();
    private static final Comparator<? super double[]> wrapComparator = new Comparator<double[]>() {

        @Override
        public int compare(double[] o1, double[] o2) {
            if (o1[0] < o2[0]) return -1;
            else if (o1[0] > o2[0]) return 1;
            else return 0;
        }
    };
    public static ErrorModel errorModel;
    private static List<List<int[]>> intArrLists;
    private static List<List<double[][]>> wrapArrLists;
    private static List<Solution1> solutionList;
    private static List<List<Solution1>> solutionLists;

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

    static {
        intArrLists = new ArrayList<List<int[]>>();
        wrapArrLists = new ArrayList<List<double[][]>>();
        solutionList = new ArrayList<Solution1>();
        solutionLists = new ArrayList<List<Solution1>>();
        errorModel = new PoissonErrorModel();
    }

    public static int[] borrowIntArray(int length) {
        int[] arr = null;
        if (intArrLists.size() <= length) {
            arr = new int[length];
        } else {
            List<int[]> intArrList = intArrLists.get(length);
            int size = intArrList.size();
            if (size == 0) {
                arr = new int[length];
            } else {
                arr = intArrList.remove(size - 1);
                Arrays.fill(arr, 0);
            }
        }
        return arr;
    }

    public static void returnIntArray(int[] arr) {
        while (arr.length >= intArrLists.size()) {
            intArrLists.add(new ArrayList<int[]>());
        }

        intArrLists.get(arr.length).add(arr);
    }

    private static double[][] borrowWrappArray(int length) {
        double[][] arr = null;
        if (wrapArrLists.size() <= length) {
            arr = new double[length][2];
        } else {
            List<double[][]> wrapArrList = wrapArrLists.get(length);
            int size = wrapArrList.size();
            if (size == 0) {
                arr = new double[length][2];
            } else {
                arr = wrapArrList.remove(size - 1);
//				Arrays.fill(arr, 0);
            }
        }
        return arr;
    }

    private static void returnWrapArray(double[][] arr) {
        while (arr.length >= wrapArrLists.size()) {
            wrapArrLists.add(new ArrayList<double[][]>());
        }

        wrapArrLists.get(arr.length).add(arr);
    }


    public static Solution1 borrowSolution(int[] counts, Signature s, int r, double weight) {
        Solution1 solution = borrowSolution(counts.length, r);
        fill(solution.counts, counts, solution.counts.length);
        solution.s.setTo(s);
        solution.weight = weight;
        return solution;
    }

    public static Solution1 borrowSolution(int[] counts, int[] s, int r, double weight) {
        Solution1 solution = borrowSolution(counts.length, r);
        fill(solution.counts, counts, solution.counts.length);
        solution.s.setTo(s);
        solution.weight = weight;
        return solution;
    }


    public static Solution1 borrowSolution(int k, int r) {
        Solution1 solution;
        if (solutionList.isEmpty()) {
            solution = new Solution1();
        } else {
            solution = solutionList.remove(solutionList.size() - 1);
        }
        solution.counts = borrowIntArray(k);
        solution.weight = 0;
        return solution;
    }

    public static void returnSolution(Solution1 solution) {

        returnIntArray(solution.counts);
//		TODO: add pooling mechanisim for trove lists
//		returnIntArray(solution.s);
        solution.counts = null;
//		solution.s = null;
        solution.weight = 0;
        solution.totalCount = 0;
        solutionList.add(solution);
    }

    public static List<Solution1> borrowSolutionList() {
        if (solutionLists.isEmpty()) {
            return new ArrayList<Solution1>();
        } else return solutionLists.remove(solutionLists.size() - 1);
    }

    public static void returnSolutionList(List<Solution1> list) {
        for (Solution1 solution : list) {
            Env.returnSolution(solution);
        }
        list.clear();
        solutionLists.add(list);
    }


    public static void fill(int[] toFill, int[] original, int last) {
        for (int i = 0; i < last; ++i) {
            toFill[i] = original[i];
        }
    }

    public static void fillSuffix(int[] toFill, int[] original) {
        for (int i = original.length - 1; i >= 0; --i) {
            toFill[i + 1] = original[i];
        }
    }

    public static final int dig(int m) {
        int d = 0;
        while (m % 2 == 0) {
            m /= 2;
            ++d;
        }
        return d;
    }

    public static int getMaxR(int[] counts) {
        return getMaxR(counts, counts.length - 1);
    }

    public static int getMaxR(Weights weights) {
        return getMaxR(weights, weights.length() - 1);
    }


    public static int getMaxR(int[] counts, int lastCount) {
        int maxValue = 1;
        for (int k = lastCount; k >= 0; --k) {
            if (maxValue < counts[k]) {
                maxValue = counts[k];
            }
        }
        return (int) Math.ceil(Math.log(maxValue + 1) / Math.log(2)) + 3;
    }

    public static int getMaxR(Weights weights, int lastCount) {
        int maxValue = 1;
        for (int k = lastCount; k >= 0; --k) {
            int maxCount = weights.getMaxCount(k);
            if (maxValue < maxCount) {
                maxValue = maxCount;
            }
        }
        return (int) Math.ceil(Math.log(maxValue + 1) / Math.log(2)) + 3;
    }

    public static void rank(List<Double> toRank, int[] ranks) {
        double[][] wrapped = borrowWrappArray(toRank.size());
        for (int i = toRank.size() - 1; i >= 0; --i) {
            wrapped[i][0] = toRank.get(i);
            wrapped[i][1] = i;
        }

        Arrays.sort(wrapped, wrapComparator);

        for (int i = toRank.size() - 1; i >= 0; --i) {
            ranks[i] = (int) (wrapped[i][1] + 0.1);
        }

        returnWrapArray(wrapped);
    }

    public static BufferedReader getBufferedReader(String path)
            throws IOException {
        if (path.endsWith(".gz")) {
            return new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(path)))));
        } else return new BufferedReader(new FileReader(path));
    }

    public static BufferedWriter getBufferedWriter(String path)
            throws IOException {
        if (path.endsWith(".gz")) {
            return new BufferedWriter(new OutputStreamWriter((new GZIPOutputStream(new FileOutputStream(path)))));
        } else return new BufferedWriter(new FileWriter(path));
    }

    public static int max(int... values) {
        int maxColumn = values[0];
        for (int i = 1; i < values.length; ++i) {
            maxColumn = Math.max(maxColumn, values[i]);
        }
        return maxColumn;
    }

    public static Pattern getPattern(int maxColumn, String sepPtrnStr, String entPtrnStr) {
        String ptternStr = entPtrnStr;
        for (int i = 1; i < maxColumn; ++i) {
            ptternStr += sepPtrnStr + entPtrnStr;
        }
        Pattern p = Pattern.compile(ptternStr);
        return p;
    }

    public static Map<String, Integer> columnTitles2Ixs(String titleLine, String separator) {
        Map<String, Integer> map = new HashMap<String, Integer>();
        String[] titles = titleLine.split(separator);
        for (int i = 0; i < titles.length; ++i) {
            map.put(titles[i], i + 1);
        }
        return map;
    }

    public static int[] count(String str) {
        int maxChar = -1;
        for (int i = str.length() - 1; i >= 0; --i) {
            maxChar = Math.max(maxChar, Character.toUpperCase(str.charAt(i)) - 'A');
        }

        int[] counts = new int[maxChar + 1];
        for (int i = str.length() - 1; i >= 0; --i) {
            ++counts[Character.toUpperCase(str.charAt(i)) - 'A'];
        }

        return counts;
    }


}
