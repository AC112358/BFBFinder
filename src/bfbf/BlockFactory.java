package bfbf;

import bfbf.Block;
import bfbf.Palindrome;

public class BlockFactory extends PalindromeFactory<Block> {

	@Override
	public Block make(Palindrome center, Palindrome wrap) {
		return new Block(center);
	}

}
