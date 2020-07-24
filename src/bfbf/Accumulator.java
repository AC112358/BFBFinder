package bfbf;

public abstract class Accumulator<In, Out> {

	public final String inValueTitle;
	public final String outValueTitle;
	
	protected Accumulator(String inValueTitle, String outValueTitle){
		this.inValueTitle = inValueTitle;
		this.outValueTitle = outValueTitle;
	}

	public abstract void accumulate(In input);
	public abstract Out accumulatedValue();
	public abstract void clear();
}
