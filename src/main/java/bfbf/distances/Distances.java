package bfbf.distances;
import bfbf.weights.Weights;


public class Distances {
    public Distances(){

    }
    public boolean compareCombinedWeights(double countWeight1, double fbrWeight1,
                                         double countWeight2, double fbrWeight2){
        return combineWeights(countWeight1, fbrWeight1) >= combineWeights(countWeight2, fbrWeight2);
    }

    public double combineWeights(double countWeight, double fbrWeight){

        return 0.5*countWeight + 0.5*fbrWeight;
    }

    public boolean firstPairDominates(double countWeight1, double fbrWeight1,
                                 double countWeight2, double fbrWeight2){
        return (firstPairStrictlyDominates(countWeight1, fbrWeight1, countWeight2, fbrWeight2) ||
                countWeight1 == countWeight2 && fbrWeight1 == fbrWeight2);
    }

    public boolean firstPairStrictlyDominates(double countWeight1, double fbrWeight1,
                double countWeight2, double fbrWeight2) {
        /*return (countWeight1 >= countWeight2 && fbrWeight1 > fbrWeight2 ||
                countWeight1 > countWeight2 && fbrWeight1 >= fbrWeight2 );*/
        return compareCombinedWeights(countWeight1, fbrWeight1, countWeight2, fbrWeight2);
    }
    public String toString(){
        return "Optimize average";
    }

}
