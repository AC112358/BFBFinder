package bfbf.weights;

public class PolynomialCountErrorModel extends ErrorModel{

    public String toString() {
        return "Low error for high values (1 - quadratic)";
    }

    @Override
    public double weight(int trueCount, int estimatedCount) {
        int diff = Math.abs(estimatedCount - trueCount);
        double threshold = -0.007 * estimatedCount * estimatedCount + 0.4 * estimatedCount - 0.4;
        //System.out.println(estimatedCount + "," + diff + ": " + (1 - (diff * diff + 0.0)/(threshold * threshold)));
        return Math.min(1, Math.max(0, 1 - 0.5 * (diff * diff + 0.0)/(threshold * threshold)));
    }
}
