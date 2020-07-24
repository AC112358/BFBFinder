package cgh;

public class EqualValues extends Accumulator<String, String> {
	
	String value = null;

	protected EqualValues(String inValueTitle, String outValueTitle) {
		super(inValueTitle, outValueTitle);
	}

	@Override
	public void accumulate(String input) {
		if (value != null){
			if (!value.equals(input)){
				throw new RuntimeException("Unequal values: expecting " + value +", obtained " + input); 
			}
		}
		else{
			value = input;
		}
	}

	@Override
	public String accumulatedValue() {
		return value;
	}

	@Override
	public void clear() {
		value = null;
	}

}
