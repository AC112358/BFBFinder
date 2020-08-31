package bfbf.weights;

import bfbf.weights.ErrorModel;

public class LowerBoundRelativeErrorModel extends ErrorModel {
    /**
     * An error model that strongly punishes deviations lower than the input counts
     * and not as much above
     *
     * @author zakov
     */

        public String toString() {
            return "Lower bound relative error";
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
                return 0.25*((double) trueCount) / estimatedCount;
            }else{
                //System.out.println(0.25*((double) estimatedCount) / trueCount);
                return ((double) estimatedCount) / trueCount;
            }
        }

//	@
//	}

    }
