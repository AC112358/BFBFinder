package cgh;

public class SumOfSquarePowers extends Accumulator<String, Double> {

	private final double base;
	private double sum;
	
	public SumOfSquarePowers(String inValueTitle, String outValueTitle, double base){
		super(inValueTitle, outValueTitle);
		this.base = base;
		sum = 0;
	}

	@Override
	public void accumulate(String input) {
		double signal = Math.pow(base, Double.parseDouble(input));
		sum += signal*signal;
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
