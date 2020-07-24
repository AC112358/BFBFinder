import java.util.List;

import bfbf.Signature;


public class SigTmp {

	private static final double factor = 1/Math.log(2);

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		System.out.println(Math.sqrt(2));
		int maxN = getMaxN();
		System.out.println(maxN);
		double logMaxN = Math.log(maxN) * factor;
		System.out.println(logMaxN);
		System.out.println(logMaxN * logMaxN) ;
		
		
		List<Long> ls = Signature.numOfSignatureLs(7000);
		double a = (Math.sqrt(5) - 1)/2;//0.62975;//  0.3575;
		double b = -3.0147001543285829594;//-1.44275;

//		a = .4;
//		b = -3;
		
		for (int i=0; i<100; ++i){
			System.out.println(i + ", " + ls.get(i) + ", " + g(a, b, 0, i) + ", " + ls.get(i)/g(a, b, 0, i));
		}
		
		for (int i=ls.size()-100; i<ls.size(); ++i){
			System.out.println(i + ", " + ls.get(i) + ", " + g(a, b, 0, i) + ", " + ls.get(i)/g(a, b, 0, i));
		}

		System.exit(0);

		double[] v = new double[8];
		int size = ls.size();
		for (int i=0; i<8; ++i){
			v[i] = ls.get(size-i-1) / f(a, b, size-i-1);
		}
		
		for (int j=0; j<2; ++j){
			for (int i=0; i<8; i+=2){
				System.out.print(v[7-i-j] + ", ");
			}
			System.out.println();
			System.out.println((v[7-j-2] - v[7-j]) + ", " + (v[7-j-6] - v[7-j-4]) +
					", " + (-v[7-j-2] + v[7-j] + v[7-j-6] - v[7-j-4]));
			System.out.println();
		}

		double v4 = ls.get(size-4) / f(a, b, size-4);
		double v2 = ls.get(size-2) / f(a, b, size-2);
		System.out.println(v4 + ", " + v2 + ", " + (v2-v4));
		
		double v3 = ls.get(size-3) / f(a, b, size-3);
		double v1 = ls.get(size-1) / f(a, b, size-1);
		System.out.println(v3 + ", " + v1 + ", " + (v1-v3));

		printLogRatios();
		//		for (int i=size-3; i<ls.size(); i+=2){
//			System.out.print(ls.get(i) / f(a, b, i) + ",");
//		}
	}
	
	private static void printLogRatios() {
		List<Long> ls = Signature.numOfSignatureLs(7336);
		for (int i=0; 2*i <= 7336; ++i){
			Long ai = ls.get(i);
			Long a2i = ls.get(2*i);
			System.out.println(i + "\t" + ai + "\t" + a2i + "\t" 
			+ Math.log(i) + "\t" + Math.log(((double) a2i)/ai)*factor);
		}
		
	}
	
	private static int getMaxN() {
		int max = 10000;
		List<Long> ls = Signature.numOfSignatureLs(max);
		int i=1;
		for (; i<=max && ((double) ls.get(i/2)) + ls.get(i-1) < Long.MAX_VALUE; ++i);
		System.out.println(ls.get(i-1) + ", " + Long.MAX_VALUE + ", " + Math.pow(2, 63) 
				+ ", " + Math.pow(2, 63)/ls.get(i-1));
		return i-1;
	}

	public static double f(double a, double b, int n){
		double logn = Math.log(n) * factor;
		return Math.pow(2, a*logn*logn + b);
	}

	public static double g(double a, double b, double c, int n){
		double logn = Math.log(n) * factor;
		return Math.pow(2, a*logn*logn + b*logn + c);
	}
}
