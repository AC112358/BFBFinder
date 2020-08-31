package main.java.bfbf.weights;

import bfbf.weights.ErrorModel;

public class LowerBoundExactErrorModel extends ErrorModel {
    public String toString() {
        return "Lower bound only";
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
        }else{
            return estimatedCount/(trueCount+0.0);
        }
    }
}
