package bfbf.weights;

public class LowerBoundPoissonErrorModel extends ErrorModel{
    public String toString() {
        return "Lower bound Poisson error";
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
        //System.out.println(trueCount + " " + estimatedCount);
        if (trueCount < estimatedCount) {
            //System.out.println(((double) trueCount) / estimatedCount);
            return 0;
        }else{
            //System.out.println(0.25*((double) estimatedCount) / trueCount);
            PoissonErrorModel pm = new PoissonErrorModel();
            return pm.weight(trueCount, estimatedCount);
        }
    }
}
