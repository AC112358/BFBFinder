package cgh;

public class SumOfPowers extends Accumulator<String, Double> {

	private final double base;
	private double sum;
	
	public SumOfPowers(String inValueTitle, String outValueTitle, double base){
		super(inValueTitle, outValueTitle);
		this.base = base;
		sum = 0;
	}

	@Override
	public void accumulate(String input) {
		sum += Math.pow(base, Double.parseDouble(input));
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
