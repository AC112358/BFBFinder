package cgh;

public class Sum extends Accumulator<String, Double> {

    private double sum;

    public Sum(String inValueTitle, String outValueTitle) {
        super(inValueTitle, outValueTitle);
        sum = 0;
    }

    @Override
    public void accumulate(String input) {
        sum += Double.parseDouble(input);
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
