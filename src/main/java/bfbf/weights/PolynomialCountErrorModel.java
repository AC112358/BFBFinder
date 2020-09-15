package bfbf.weights;

public class PolynomialCountModel extends ErrorModel{

    public String toString() {
        return "Low error for high values";
    }

    @Override
    public double weight(int trueCount, int estimatedCount) {
        
    }
}
