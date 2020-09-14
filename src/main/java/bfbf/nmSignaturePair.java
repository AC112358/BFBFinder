package bfbf;

import java.util.Collections;

import bfbf.Signature;

public class nmSignaturePair implements Comparable<nmSignaturePair>{
    protected Signature sig;
    protected int nmSum;

    public nmSignaturePair(int nmSum, Signature sig){
        this.nmSum = nmSum;
        this.sig = sig;
    }

    public int compareTo(nmSignaturePair other) {
        if (nmSum != other.nmSum){
            return nmSum - other.nmSum;
        }
        if (sig == null || other.sig == null){
            return 0;
        }
        return sig.compareTo(other.sig);
    }

    @Override
    public String toString(){
        return "(" + nmSum + " ; " + sig.toString() + ")";
    }
}
