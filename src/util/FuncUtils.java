package util;

import java.util.Random;

/**
 * 
 * @author: Yang Qian(HeFei University of Technology)
 */
public class FuncUtils {
	/**
	 * Sample a value from a double array
	 * 
	 * @param probs 
	 * @return
	 */
	public static int rouletteGambling(double[] prob){
		int topic = 0;
		for (int k = 1; k < prob.length; k++) {
			prob[k] += prob[k - 1];
		}
		double u = Math.random() * prob[prob.length - 1];
		for (int t = 0; t < prob.length; t++) {
			if (u < prob[t]) {
				topic = t;
				break;
			}
		}
		return topic;
	}
	/**
	 * Sample a value from a double array
	 * 
	 * @param probs 
	 * @return
	 */
	public static int rouletteGambling(double[][] prob){
		int K = prob.length;
		int A = prob[0].length;
		double[] pr_sum = new double[K * A];
		for (int k = 0; k < K; k++) {
			for (int a = 0; a < A; a++) {
				pr_sum[k  + a*K] = prob[k][a];
			}
		}
		int idx = rouletteGambling(pr_sum);
		return idx;
	}
	/**
	 * get the index of max value in an array
	 * 
	 * @param array
	 * @return the index of max value
	 */
	public static int maxValueIndex(double[] array) {
		double max = array[0];
		int maxVIndex = 0;
		for (int i = 1; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
				maxVIndex = i;
			}
		}
		return maxVIndex;
	}

	public static double[] getGaussianSample(int K, double mean, double deviation) {
		Random r = new Random();
		double[] sample = new double[K];
		for(int k = 0; k < K; k ++) {
			sample[k] = r.nextGaussian() * Math.sqrt(deviation) + mean;
		}
		return sample;
	}
}
