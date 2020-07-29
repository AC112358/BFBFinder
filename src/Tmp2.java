public class Tmp2 {

    public static void main(String[] args) {
        int k = 1025;
        long[] series = new long[k];
        series[0] = 1;
        series[1] = 1;

        for (int n = 2; n < k; ++n) {
            series[n] = series[n - 1] + series[n / 2];
            System.out.println(n + ", " + series[n] + ", "
                    + an(n) + ", " + series[n] / an(n) + ", " + (an(n) >= series[n]));
        }
    }

    static double an(int n) {
        double logn = Math.floor(Math.log(n) / Math.log(2));
//		double exp = logn*logn/2;
        double exp = logn * logn / 2 + 3 / 4 * logn;
        return Math.pow(2, exp);
    }

}
