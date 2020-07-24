package bfbf;

import bfbf.Block;
import bfbf.CompositPalindrome;
import bfbf.Palindrome;

public class CompositeFactory extends PalindromeFactory<CompositPalindrome> {

	@Override
	public CompositPalindrome make(Palindrome center, Palindrome wrap) {
		return new CompositPalindrome((ConvexPalindrome) center, (Block) wrap);
	}

}
