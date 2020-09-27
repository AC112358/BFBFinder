package bfbf.weights;

public class SlightLowerBoundErrorModel extends ErrorModel{
    int upperBound = 7;
    public String toString() {
        return "Lower bound only, diff <= " + upperBound;
    }

        /*@Override
        public int minRealValue(int observedValue, double maxCurrError) {

            return observedValue;
        }

        @Override
        public int maxRealValue(int observedValue, double maxCurrError) {

            return observedValue;
        }*/


    @Override
    public double weight(int trueCount, int estimatedCount) {
        if (trueCount < estimatedCount) {
            return 0;
        }
        if (trueCount - estimatedCount <= upperBound){
            return 1;
        }
        return 0;
    }
}
