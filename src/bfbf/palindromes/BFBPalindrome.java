package bfbf.palindromes;

public abstract class BFBPalindrome extends Palindrome {

	protected BFBPalindrome(Palindrome center, Palindrome wrap) {
		super(center, wrap);
	}
	
	public Block wrap(){
		return Block.make(this);
	}


}
