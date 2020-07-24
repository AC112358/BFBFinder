package bfbf.palindromes;

public abstract class PalindromeFactory<P extends Palindrome> {
	
	public abstract P make(Palindrome center, Palindrome wrap);

}
