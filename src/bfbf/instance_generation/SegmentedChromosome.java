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

package bfbf.instance_generation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A class representing segmented chromosomal arms.
 *
 * @author Shay Zakov
 */
public class SegmentedChromosome {

    private final int ploidy;
    private List<Long> segmentLengths;
    private List<Integer>[] ploidSegments;
    private long[] ploidLengths;
    private int breakedSegmentIx;
    private long firstSubSegmentLength;
    private long originalLength;

    @SuppressWarnings("unchecked")
    public SegmentedChromosome(long length, int ploidy) {
        originalLength = length;
        this.ploidy = ploidy;
        segmentLengths = new ArrayList<Long>();
        ploidSegments = new List[ploidy];
        ploidLengths = new long[ploidy];
        for (int i = 0; i < ploidy; ++i) {
            ploidSegments[i] = new ArrayList<Integer>();
        }
        clear();
    }

    public void clear() {
        segmentLengths.clear();
        segmentLengths.add(0l);
        segmentLengths.add(originalLength);
        for (int i = 0; i < ploidy; ++i) {
            ploidSegments[i].clear();
            ploidSegments[i].add(1);
            ploidLengths[i] = originalLength;
        }
    }

    public long getBreakpoint(int ploidIx, double relativeBreakPosition) {
        return getBreakpoint(ploidIx, relativeBreakPosition, 0);
    }

    public long getBreakpoint(int ploidIx, double relativeBreakPosition, long lengthReduction) {
        return (long) (relativeBreakPosition * (ploidLength(ploidIx) - lengthReduction));
    }

    public void setPloidBreakPoint(int ploidIx, double relativeBreakPosition) {
        setPloidBreakPoint(ploidIx, getBreakpoint(ploidIx, relativeBreakPosition));
    }

    public void setPloidBreakPoint(int ploidIx, long breakPosition) {
        long accumulatedLength = 0;
        long nextSegmentLength = segmentLength(ploidSegments[ploidIx].get(0));
        for (breakedSegmentIx = 0; accumulatedLength + nextSegmentLength <= breakPosition; ++breakedSegmentIx) {
            accumulatedLength = accumulatedLength + nextSegmentLength;
            nextSegmentLength = segmentLength(ploidSegments[ploidIx].get(breakedSegmentIx + 1));
        }
        firstSubSegmentLength = breakPosition - accumulatedLength;
    }

    public void setRefBreakPoint(long breakPosition) {
        long accumulatedLength = 0;
        long nextSegmentLength = segmentLength(1);
        for (breakedSegmentIx = 1; breakedSegmentIx < segmentLengths.size() &&
                accumulatedLength + nextSegmentLength <= breakPosition; ++breakedSegmentIx) {
            accumulatedLength = accumulatedLength + nextSegmentLength;
            nextSegmentLength = segmentLength(breakedSegmentIx + 1);
        }
        firstSubSegmentLength = breakPosition - accumulatedLength;
    }


//	private void breakRef(long breakPosition){
//		setRefBreakPoint(breakPosition);
//		breakSegment(breakedSegmentIx, firstSubSegmentLength);
//	}

    private void breakPloid(int ploidIx, long breakPosition) {
        setPloidBreakPoint(ploidIx, breakPosition);
        int segment = ploidSegments[ploidIx].get(breakedSegmentIx);
        if (segment > 0) {
            breakSegment(segment, firstSubSegmentLength);
        } else {
            breakSegment(-segment, segmentLength(segment) - firstSubSegmentLength);
        }
    }

    private void breakSegment(int breakedSegment, long prefixLength) {
        long breakedSegmentLength = segmentLength(breakedSegment);
        long suffixLength = breakedSegmentLength - prefixLength;

        if (prefixLength == 0 || suffixLength == 0) {
            return;
        }

        segmentLengths.set(breakedSegment, prefixLength);
        segmentLengths.add(breakedSegment + 1, suffixLength);

        for (int ploidIx = 0; ploidIx < ploidy; ++ploidIx) {
            for (int i = 0; i < ploidSegments[ploidIx].size(); ++i) {
                int currSegment = ploidSegments[ploidIx].get(i);
                int sign = (int) Math.signum(currSegment);

                if (sign * currSegment > breakedSegment) {
                    ploidSegments[ploidIx].set(i, currSegment + sign);
                } else if (sign * currSegment == breakedSegment) {
                    ploidSegments[ploidIx].add(i + (1 + sign) / 2, currSegment + sign);
                    ++i;
                }
            }
        }
    }

    public void delete(int ploidIx, double relBreakpoint1, double relBreakpoint2) {
        delete(ploidIx, getBreakpoint(ploidIx, relBreakpoint1), getBreakpoint(ploidIx, relBreakpoint2));
    }

    public void delete(int ploidIx, long breakpoint, long segmentLength) {
        breakPloid(ploidIx, breakpoint);
        breakPloid(ploidIx, breakpoint + segmentLength);
        setPloidBreakPoint(ploidIx, breakpoint);
        int start = breakedSegmentIx;
        setPloidBreakPoint(ploidIx, breakpoint + segmentLength);
        int end = breakedSegmentIx;

        for (int i = end - 1; i >= start; --i) {
            ploidSegments[ploidIx].remove(i);
        }
        setTotalLength(ploidIx, ploidLength(ploidIx) - segmentLength);
    }

    public void tandemDuplicate(int ploidIx, double relBreakpoint1, double relBreakpoint2, boolean isInverted) {
        tandemDuplicate(ploidIx, getBreakpoint(ploidIx, relBreakpoint1), getBreakpoint(ploidIx, relBreakpoint2), isInverted);
    }

    public void tandemDuplicate(int ploidIx, long breakPosition, long segmentLength, boolean isInverted) {
        duplicate(ploidIx, breakPosition, segmentLength, breakPosition + segmentLength, isInverted);
    }

    public void duplicate(int ploidIx, double relBreakpoint1, double relBreakpoint2, double relInsertion, boolean isInverted) {
        long firstBreakPosition = getBreakpoint(ploidIx, relBreakpoint1);
        long secondBreakPosition = getBreakpoint(ploidIx, relBreakpoint2);
        if (secondBreakPosition < firstBreakPosition) {
            long tmp = firstBreakPosition;
            firstBreakPosition = secondBreakPosition;
            secondBreakPosition = tmp;
        }

        long insertionPosition = (long) ((ploidLength(ploidIx) - secondBreakPosition + firstBreakPosition) * relInsertion);
        if (insertionPosition >= firstBreakPosition) {
            insertionPosition += secondBreakPosition - firstBreakPosition;
        }
        duplicate(ploidIx, firstBreakPosition, secondBreakPosition, insertionPosition, isInverted);
    }

    public void duplicate(int ploidIx, long breakPosition, long segmentLength, long insertionPosition, boolean isInverted) {
        breakPloid(ploidIx, breakPosition);
        breakPloid(ploidIx, breakPosition + segmentLength);
        breakPloid(ploidIx, insertionPosition);
        setPloidBreakPoint(ploidIx, breakPosition);
        int start = breakedSegmentIx;
        setPloidBreakPoint(ploidIx, breakPosition + segmentLength);
        int end = breakedSegmentIx;
        setPloidBreakPoint(ploidIx, insertionPosition);
        int insertionIx = breakedSegmentIx;
        setTotalLength(ploidIx, ploidLength(ploidIx) + segmentLength);
        insert(ploidIx, start, end, insertionIx, isInverted);
    }

    public void translocate(int ploidIx, double relBreakpoint1, double relBreakpoint2, double relInsertion, boolean isInverted) {
        long firstBreakPosition = getBreakpoint(ploidIx, relBreakpoint1);
        long secondBreakPosition = getBreakpoint(ploidIx, relBreakpoint2);
        if (secondBreakPosition < firstBreakPosition) {
            long tmp = firstBreakPosition;
            firstBreakPosition = secondBreakPosition;
            secondBreakPosition = tmp;
        }

        long insertionPosition = (long) ((ploidLength(ploidIx) - secondBreakPosition + firstBreakPosition) * relInsertion);
        if (insertionPosition >= firstBreakPosition) {
            insertionPosition += secondBreakPosition - firstBreakPosition;
        }
        translocate(ploidIx, firstBreakPosition, secondBreakPosition, insertionPosition, isInverted);
    }

    public void translocate(int ploidIx, long breakpoint, long segmentLength, long insertionPosition, boolean isInverted) {
        breakPloid(ploidIx, breakpoint);
        breakPloid(ploidIx, breakpoint + segmentLength);
        breakPloid(ploidIx, insertionPosition);
        setPloidBreakPoint(ploidIx, breakpoint);
        int start = breakedSegmentIx;
        setPloidBreakPoint(ploidIx, breakpoint + segmentLength);
        int end = breakedSegmentIx;
        setPloidBreakPoint(ploidIx, insertionPosition);
        int insertionIx = breakedSegmentIx;
        insert(ploidIx, start, end, insertionIx, isInverted);
        if (insertionIx <= start) {
            int shift = end - start;
            start += shift;
            end += shift;
        }
        for (int i = end - 1; i >= start; --i) {
            ploidSegments[ploidIx].remove(i);
        }
    }

//	public void invert(int ploidIx, double relBreakpoint1, double relBreakpoint2) {
//		invert(ploidIx, getBreakpoint(ploidIx, relBreakpoint1), getBreakpoint(ploidIx, relBreakpoint2));
//	}
//

    public void invert(int ploidIx, long breakpoint, long segmentLength) {
        breakPloid(ploidIx, breakpoint);
        breakPloid(ploidIx, breakpoint + segmentLength);
        setPloidBreakPoint(ploidIx, breakpoint);
        int start = breakedSegmentIx;
        setPloidBreakPoint(ploidIx, breakpoint + segmentLength);
        int end = breakedSegmentIx;
        invert(ploidSegments[ploidIx], start, end);
    }


    private void insert(int ploidIx, int start, int end,
                        int insertionIx, boolean isInverted) {
        List<Integer> insert = new ArrayList<Integer>(ploidSegments[ploidIx].subList(start, end));
        if (isInverted) {
            invert(insert, 0, end - start);
        }
        ploidSegments[ploidIx].addAll(insertionIx, insert);
    }

    private void invert(List<Integer> ls, int start, int end) {
        for (int i = (end - start - 1) / 2; i >= 0; --i) {
            int tmp = ls.get(start + i);
            ls.set(start + i, -ls.get(end - i - 1));
            ls.set(end - i - 1, -tmp);
        }
    }

//	public void bfbCycle(int ploidIx, double relativeBreakPosition){
//		bfbCycle(ploidIx, getBreakpoint(ploidIx, relativeBreakPosition));
//	}


    public void bfbCycle(int ploidIx, long suffixLength) {
        long ploidLength = ploidLength(ploidIx);
        long breakPosition = ploidLength - suffixLength;
        breakPloid(ploidIx, breakPosition);
        setPloidBreakPoint(ploidIx, breakPosition);
        setTotalLength(ploidIx, ploidLength + suffixLength);
        for (int i = ploidSegments[ploidIx].size() - 1; i >= breakedSegmentIx; --i) {
            ploidSegments[ploidIx].add(-ploidSegments[ploidIx].get(i));
        }
    }

    private void setTotalLength(int ploidIx, long newTotalLength) {
        ploidLengths[ploidIx] = newTotalLength;
    }

    private long segmentLength(int segment) {
        return segmentLengths.get(Math.abs(segment));
    }


    public int[] getCounts() {
        int[] counts = new int[segmentLengths.size()];
        for (int ploidIx = 0; ploidIx < ploidy; ++ploidIx) {
            for (int i = 0; i < ploidSegments[ploidIx].size(); ++i) {
                ++counts[Math.abs(ploidSegments[ploidIx].get(i))];
            }
        }
        return counts;
    }

    public boolean unifyConsecutiveDeletedSegments() {
        int[] counts = getCounts();

        List<Integer> toRemove = null;

        for (int i = 1; i < counts.length; ++i) {
            if (counts[i] == 0) {
                int j = i + 1;
                while (j < counts.length && counts[j] == 0) {
                    if (toRemove == null) {
                        toRemove = new ArrayList<Integer>();
                    }
                    toRemove.add(j);
                    segmentLengths.set(i, segmentLength(i) + segmentLength(j));
                    ++j;
                }
                i = j;
            }
        }

        if (toRemove != null) {
            removeSegments(toRemove);
            for (int i = toRemove.size() - 2; i >= 0; --i) {
                segmentLengths.remove((int) toRemove.get(i));
            }
        }
        return toRemove != null;
    }

    public String toString() {
        String res = "";
        for (int ploidIx = 0; ploidIx < ploidy; ++ploidIx) {
            res += "chromosome " + ploidIx + ": " + ploidSegments[ploidIx].toString() + "\n";
        }
        res += "Counts: " + Arrays.toString(getCounts())
                + "\nSegment lengths: " + segmentLengths
                + "\nDiscordant breakpoints: " + discordantBreakpoints()
                + ", foldback breakpoints: " + foldbackBreakpoints()
//												+ "\nTotal weight error: " + (calcAllPloidLengths() - allPloidLength()) + 
//												", segment weight error:" + (originalLength - calcTotalSegmentLength())
        ;
        return res.toString();
    }

    public int discordantBreakpoints() {
        int discordantBreakpoints = 0;

        for (int ploidIx = 0; ploidIx < ploidy; ++ploidIx) {
            int currSegment = ploidSegments[ploidIx].get(0), nextSegment;
            for (int i = 1; i < ploidSegments[ploidIx].size(); ++i) {
                nextSegment = ploidSegments[ploidIx].get(i);
                if (nextSegment != currSegment + 1) {
                    ++discordantBreakpoints;
                }
                currSegment = nextSegment;
            }
        }
        return discordantBreakpoints;
    }

    public int foldbackBreakpoints() {
        int foldbackBreakpoints = 0;
        for (int ploidIx = 0; ploidIx < ploidy; ++ploidIx) {
            int currSegment = ploidSegments[ploidIx].get(0), nextSegment;
            for (int i = 1; i < ploidSegments[ploidIx].size(); ++i) {
                nextSegment = ploidSegments[ploidIx].get(i);
                if (nextSegment == -currSegment) {
                    ++foldbackBreakpoints;
                }
                currSegment = nextSegment;
            }
        }
        return foldbackBreakpoints;
    }

    public long ploidLength(int ploidIx) {
        return ploidLengths[ploidIx];
    }

    @Deprecated
    public boolean mergeShortSegments(long minSegmentLength) {
        if (minSegmentLength == 0) return false;
        List<Integer> toRemove = null;

        for (int i = 1; i < segmentLengths.size() - 1; ++i) {
            if (segmentLengths.get(i) < minSegmentLength) {
                if (toRemove == null) {
                    toRemove = new ArrayList<Integer>();
                }
                toRemove.add(i);
                segmentLengths.set(i + 1, segmentLengths.get(i) + segmentLengths.get(i + 1));
            }
        }

        if (segmentLengths.get(segmentLengths.size() - 1) < minSegmentLength) {
            if (toRemove == null) {
                toRemove = new ArrayList<Integer>();
            }
            toRemove.add(segmentLengths.size() - 1);
            int i = toRemove.size() - 2;
            while (toRemove.get(i) + 1 == toRemove.get(i + 1)) {
                --i;
            }
            segmentLengths.set(toRemove.get(i + 1) - 1, segmentLengths.get(toRemove.get(i + 1) - 1) +
                    segmentLengths.get(segmentLengths.size() - 1));
        }

        if (toRemove != null) {
            removeSegments(toRemove);
            return true;
        } else return false;
    }

    private void removeSegments(List<Integer> toRemove) {
        toRemove.add(segmentLengths.size());

        int[] fixedSegments = new int[segmentLengths.size()];
        for (int i = 0, segment = 1; i < toRemove.size(); ++i) {
            int r = toRemove.get(i);
            for (; segment < r; ++segment) {
                fixedSegments[segment] = segment - i;
            }
        }

        List<Integer> list = new ArrayList<Integer>();
        for (int ploidIx = 0; ploidIx < ploidy; ++ploidIx) {
            for (int segment : ploidSegments[ploidIx]) {
                int sign = (int) Math.signum(segment);
                if (!toRemove.contains(sign * segment)) {
                    list.add((int) (fixedSegments[sign * segment] * sign));
                }
            }
            ploidSegments[ploidIx].clear();
            ploidSegments[ploidIx].addAll(list);
            list.clear();
        }
    }

    public List<Integer>[] getPloidSegments() {
        return ploidSegments;
    }

    public void keepRefSegmentation() {
        for (int i = 0; i < ploidy; ++i) {
            ploidSegments[i].clear();
            for (int j = 1; j < segmentLengths.size(); ++j) {
                ploidSegments[i].add(j);
            }
            ploidLengths[i] = originalLength;
        }
    }


    //	public List<Integer> segmentList() {
    //		return ploidSegments;
    //	}
    //

    //	public void mergeBreakpoints(SegmentedChromosome other) {
    //		int position = 0;
    //		for (long length : other.segmentLengths){
    //			position += length;
    //			if (position >= originalLength) break;
    //			breakRef(position);
    //		}
    //	}

}
