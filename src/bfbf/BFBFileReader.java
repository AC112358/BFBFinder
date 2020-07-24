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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import bfbf.Env;

/**
 * A reader for count vector files.
 * 
 * @author Shay Zakov
 *
 */
public class BFBFileReader {
	
	public static void readData(List<String>[] metaData, List<int[]>[] counts,
			String path) throws IOException {
		BufferedReader bf = new BufferedReader(new FileReader(path));
		String line = bf.readLine();
		String[] split, countArr;
		int[] count;
	
		while (line != null){
			split = line.split("\\s");
			assert split.length == 3 : "Invalid line: " + line;
	
			if (split[2].contains("C,")){
				metaData[0].add(split[0]);
				metaData[1].add(split[1]);
	
				split = split[2].split("C,");
	
				countArr = split[0].split(",");
				count = new int[countArr.length];
				for (int i=0; i<countArr.length; ++i){
					count[i] = Integer.parseInt(countArr[countArr.length - i - 1]);// - 1;
				}
				counts[0].add(standartize(count));
	
				countArr = split[1].split(",");
				count = new int[countArr.length];
				for (int i=0; i<countArr.length; ++i){
					count[i] = Integer.parseInt(countArr[i]);
				}
				counts[1].add(count);
			}
	
			line = bf.readLine();
		}
		bf.close();
	}

	public static List<ChromosomeArm> readData(String path) throws IOException {
		List<ChromosomeArm> arms = new ArrayList<ChromosomeArm>();

		BufferedReader bf = new BufferedReader(new FileReader(path));
		String line;
		do{
			line = bf.readLine();
		} while (line == "");
		
		String[] split, armsCounts, countArr;
		int[] count;
	
		while (line != null){
			
			split = line.split("\\s+");
			assert split.length == 3 : "Invalid line: " + line;
			
			armsCounts = split[2].replaceAll("\"", "").split("C,");
			
			for (int j=0; j<armsCounts.length; ++j){
				countArr = armsCounts[j].split("[\",]");
				count = new int[countArr.length];
				int factor = 2*j - 1;
				int shift = (1-j)*(countArr.length-1);
				for (int i=0; i<countArr.length; ++i){
					count[factor * i + shift] =  Integer.parseInt(countArr[i]);
				}
				arms.add(new ChromosomeArm(split[0], Integer.parseInt(Env.chrStr2Int.get(split[1])), j == 0, standartize(count)));
			}

			do{
				line = bf.readLine();
				if (line != null){
					line = line.trim();
				}
			} while (line != null && line.length() == 0);
		}
		bf.close();
		
		return arms;
	}

	private static int[] standartize(int[] count) {
		
		List<Integer> ls = new ArrayList<Integer>();
		int prev = 0;
		
		for (int i=0; i<count.length; ++i){
			if (count[i] != prev){
				ls.add(count[i]);
				prev = count[i];
			}
		}
		
		if (!ls.isEmpty() && ls.get(ls.size()-1) == 0){
			ls.remove(ls.size()-1);
		}
		
		if (ls.size() != count.length){
			count = new int[ls.size()];
			for (int i=0; i<count.length; ++i){
				count[i] = ls.get(i);
			}
		}
		
		return count;
		
//		int i=0, j;
//		for (; i<count.length; ++i){
//			if (count[i] < 0){
//				return new int[0];
//			}
//		}
//		for (i=0; i<count.length && count[i] == 0; ++i);
//		for (j = count.length; j > i && count[j-1] == 0; --j);
//		if (i != 0 || j != count.length){
//			return Arrays.copyOfRange(count, i, j);
//		}
//		else return count;
	}
	
	public static List<int[]> readSimCounts(String path) throws IOException{
		List<int[]> counts = new ArrayList<int[]>();
		BufferedReader bf = new BufferedReader(new FileReader(path));
		String line;
		
		String[] split, armsCounts, countArr;
		int[] count;
	
		line = bf.readLine();
		while (line != null){
			if (line.startsWith("[0, ")){
				split = line.substring(4, line.length()-1).split(",\\s*");
				count = new int[split.length];
				for (int i=0; i<split.length; ++i){
					count[i] = Integer.parseInt(split[i]);
				}
				counts.add(count);
			}
			line = bf.readLine();
		}
		
		return counts;
	}

}
