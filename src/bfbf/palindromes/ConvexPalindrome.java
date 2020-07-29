package bfbf.palindromes;

import java.util.ArrayList;
import java.util.List;

public class ConvexPalindrome extends Palindrome {

    private static final List<Palindrome> pool = new ArrayList<Palindrome>();
    private static ConvexPalindrome cashKey = new ConvexPalindrome(null, null);
    private static PalindromeFactory<ConvexPalindrome> factory = new ConvexFactory();

    protected ConvexPalindrome(BFBPalindrome center, ConvexPalindrome wrap) {
        super(center, wrap);
    }

    public static ConvexPalindrome make(BFBPalindrome center, Palindrome wrap) {
        if (depthDiff(center, wrap) <= 0) {
            throw new IllegalArgumentException("Center depth must be " +
                    "greater than wrap depth.");
        }
        if (wrap instanceof BFBPalindrome) {
            wrap = make((BFBPalindrome) wrap, null);
        }
        return make(center, wrap, cashKey, factory);
    }

    public static int compare(ConvexPalindrome first, ConvexPalindrome second) {
        if (isEmpty(first) && isEmpty(second)) return 0;
        else if (isEmpty(first)) return -1;
        else if (isEmpty(second)) return 1;
        else {
            int res = first.center.compareTo(second.center);
            if (res == 0) {
                res = compare((ConvexPalindrome) first.wrap, (ConvexPalindrome) second.wrap);
            }
            return res;
        }
    }

    public ConvexPalindrome extend(Palindrome toAdd) {
        if (isEmpty(toAdd)) return this;
        else {
            if (wrap == null) {
                return make((BFBPalindrome) center, toAdd);
            } else {
                return make((BFBPalindrome) center,
                        ((ConvexPalindrome) wrap).extend(toAdd));
            }
        }
    }

    //	@Override
//	protected int myCompareTo(Palindrome other) {
//		int res = center.compareTo(other.center);
//		if (res == 0){
//			if (wrap == null){
//				if (other.wrap != null){
//					res = -1;
//				}
//			}
//			else if (other.wrap == null){
//				res = 1;
//			}
//			else{
//				res = wrap.compareTo(other.wrap);
//			}
//		}
//		return res;
//	}

    @Override
    protected int classCompareTo(Palindrome other) {
        if (other instanceof ConvexPalindrome) return super.classCompareTo(other);
        else if (other instanceof CompositPalindrome) return -1;
        else return 1;
    }

    public int nestingDeg() {
        if (wrap == null) {
            if (center == null || EmptyPalindrome.SINGLETON.equals(center)) {
                return 0;
            } else return 1;
        } else {
            return 1 + ((ConvexPalindrome) wrap).nestingDeg();
        }
    }

    public int centerDeg() {
        return center.centerDeg();
    }

    @Override
    public int minPartDepth() {
        if (isEmpty(wrap)) return center.depth();
        else return wrap.minPartDepth();
    }

    @Override
    protected Palindrome primary() {
        return center;
    }

    @Override
    protected Palindrome secondary() {
        return wrap;
    }

    @Override
    public List<Palindrome> pool() {
        return pool;
    }

    public static class ConvexFactory extends PalindromeFactory<ConvexPalindrome> {

        @Override
        public ConvexPalindrome make(Palindrome center, Palindrome wrap) {
            return new ConvexPalindrome((BFBPalindrome) center, (ConvexPalindrome) wrap);
        }

    }
}
