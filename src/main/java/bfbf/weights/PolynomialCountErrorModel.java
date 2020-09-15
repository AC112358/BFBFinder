package bfbf.weights;

public class PolynomialCountErrorModel extends ErrorModel{

    public String toString() {
        return "Low error for high values (1/polynomial)";
    }

    @Override
    public double weight(int trueCount, int estimatedCount) {
        int diff = Math.abs(estimatedCount - trueCount);
        return Math.max(0, 1 - 10 * Math.pow(diff, 3) / Math.pow(trueCount + 1, 3));
    }
}
