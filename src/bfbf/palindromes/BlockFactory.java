package bfbf.palindromes;

public class BlockFactory extends PalindromeFactory<Block> {

	@Override
	public Block make(Palindrome center, Palindrome wrap) {
		return new Block(center);
	}

}
