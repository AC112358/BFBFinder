package cgh;

public class SumOfSquares extends Accumulator<String, Double> {
	
	private double sum;
	
	public SumOfSquares(String inValueTitle, String outValueTitle){
		super(inValueTitle, outValueTitle);
		sum = 0;
	}

	@Override
	public void accumulate(String input) {
		double signal = Double.parseDouble(input);
		sum += signal * signal;
	}

	@Override
	public Double accumulatedValue() {
		return sum;
	}

	@Override
	public void clear() {
		sum = 0;
	}

}
