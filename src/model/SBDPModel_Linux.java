package model;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.special.Gamma;

import evaluation.EstimationUtil;
import util.FileUtil;
import util.FuncUtils;

/**
 * 
 * This model can be used to cluster short text.
 * 
 * @author: Yang Qian,Yezheng Liu,Yuanchun Jiang (HeFei University of Technology)
 */
public class SBDPModel_Linux {
	public double alpha; // Hyper-parameter alpha
	//	public double beta; // Hyper-parameter beta
	public int K; // number of clusters
	public int iterations; // number of Gibbs sampling iterations
	public int topWords; // number of most probable words for each cluster
	public Map<String, Integer> wordToIndexMap = new HashMap<String, Integer>();;  //word to index
	public List<String> indexToWordMap = new ArrayList<String>();    //index to String word 
	public int [][] docword;//word index array
	public int M; // number of documents in the corpus
	public int V; // number of words in the corpus
	//DPMM
	public int[] z; //cluster assignment for each document
	public int[] kd; //number of documents assigned to a cluster
	public int[][] nkw; //cluster-word count
	public int[] nksum; //cluster-word sum (total number of words assigned to a cluster)
	//Example: given a document of "a a b a b c d c". We have: 1 2 1 3 2 1 1 2
	public List<List<Integer>> docWordIndexCount;
	public double pi_b[]; 
	public boolean b[][];  
	public int b_sum[];  
	public double beta0 = 1E-12; 
	public double beta1 = 0.1;  
	public double s = 1.0; 
	public double r = 1.0;  
	public JDKRandomGenerator rand; 
	//output
	public int topWordsOutputNumber;
	public String outputFileDirectory; 

	public SBDPModel_Linux(String inputFile, String inputFileCode, int initClusterNumber,
			double inputAlpha, double inputBeta, int inputIterations, int inTopWords,
			String outputFileDir){
		//read data
		ArrayList<String> docLines = new ArrayList<String>();
		FileUtil.readLines(inputFile, docLines,inputFileCode);
		M = docLines.size();
		docword = new int[M][];
		docWordIndexCount = new ArrayList<List<Integer>>();
		int j = 0;
		for(String line : docLines){
			List<String> words = new ArrayList<String>();
			FileUtil.tokenizeAndLowerCase(line, words);
			docword[j] = new int[words.size()];
			HashMap<String, Integer> wordCount = new HashMap<String, Integer>();
			List<Integer> wordIndexCount = new ArrayList<Integer>();
			for(int i = 0; i < words.size(); i++){
				String word = words.get(i);
				if(!wordToIndexMap.containsKey(word)){
					int newIndex = wordToIndexMap.size();
					wordToIndexMap.put(word, newIndex);
					indexToWordMap.add(word);
					docword[j][i] = newIndex;
				} else {
					docword[j][i] = wordToIndexMap.get(word);
				}
				int times = 0;
				if (wordCount.containsKey(word)) {
					times = wordCount.get(word);
				}
				times += 1;
				wordCount.put(word, times);
				wordIndexCount.add(times);
			}
			docWordIndexCount.add(wordIndexCount);
			j++;
		}
		V = indexToWordMap.size();
		alpha = inputAlpha;
		K = initClusterNumber;  //initialize topic number
		iterations = inputIterations;
		topWordsOutputNumber = inTopWords;
		outputFileDirectory = outputFileDir;
		initialize();
	}
	/**
	 * Randomly initialize topic assignments
	 */
	public void initialize(){
		rand = new JDKRandomGenerator();
		rand.setSeed(System.currentTimeMillis());
		BetaDistribution betaDist = new BetaDistribution(rand, 1.0 , 1.0);
		b = new boolean[K][V]; 
		b_sum = new int[K]; 
		pi_b = new double[K]; 
		for (int k2 = 0; k2 != K; k2++) {

			pi_b[k2] = betaDist.sample();
			for (int v = 0; v != V; v++) {
				b[k2][v] = true;
			}
			b_sum[k2] = V;
		}

		z = new int[M];
		nkw = new int[K][V];
		nksum = new int[K];
		kd = new int[K];
		for (int d = 0; d < M; d++) {
			int cluster = (int) (Math.random() * K);
			z[d] = cluster; 
			updateCount(d,cluster,false);
		}
	}
	public void MCMCSampling(){
		for (int iter = 1; iter <= iterations; iter++) {
			System.out.println("iteration : " + iter);
			gibbsOneIteration();
			defragment();
			System.out.println("K is:\t" + K);
			sampleBinaryBMatrix();
		}
		// output the result
		System.out.println("write cluster word ..." );
		writeTopWordsWithProbability();
		System.out.println("write cluster distribution ..." );
		writeClusterDistri();
		System.out.println("write cluster for each document ..." );
		writeDocCluster();
		System.out.println("write average coherence ..." );
		writeCoherence("/home/qianyang/yelp/reviews_small");

	}
	public void gibbsOneIteration() {
		for (int d = 0; d < M; d++) {
			double[] p = new double[K + 1];
			int cluster = z[d];
			//decrease the count
			updateCount(d,cluster,true);
			if (kd[cluster] < 0) {
				defragment();
			}
			//sample the new cluster
			//			System.out.println(b.length + "\t" + b[0].length);
			for (int k = 0; k < K; k++) {
				p[k] = (kd[k])/(M - 1 + alpha );
				for (int n = 0; n < docword[d].length; n ++) {
					int y = b[k][docword[d][n]] ? 1 : 0;
					p[k] *= (nkw[k][docword[d][n]] + y*beta1 + beta0 + docWordIndexCount.get(d).get(n) - 1)
							/ (nksum[k] + V*beta0 + b_sum[k]*beta1 + n);
				}
			}
			p[K] = (alpha)/(M - 1 +  alpha);
			for (int n = 0; n < docword[d].length; n ++) {
				p[K] *= (beta1 + docWordIndexCount.get(d).get(n) - 1)
						/ (V * beta1 + n);
			}
			cluster = FuncUtils.rouletteGambling(p);
			//increase the count
			z[d] = cluster;
			if (cluster < K) { //old cluster
				updateCount(d,cluster,false);
			}else {  //new cluster
				K++;
				kd = ensureCapacity(kd,1);
				nkw = ensureCapacity(nkw, 1, 0);
				nksum = ensureCapacity(nksum, 1);
				b = ensureCapacity(b, 1, 0);
				b_sum = ensureCapacity(b_sum,1);
				pi_b = ensureCapacity(pi_b,1);
				updateCount(d,cluster,false);
			}
		}
	}
	//remove the topic without words
	public void defragment() {
		int[] kOldToKNew = new int[K];
		int newK = 0;
		for (int k = 0; k < K; k++) {
			if (kd[k] > 0) {
				kOldToKNew[k] = newK;  //
				swap(kd, newK, k);
				swap(nkw, newK, k);
				swap(nksum, newK, k);
				swap(b, newK, k);
				swap(b_sum, newK, k);
				swap(pi_b, newK, k);
				newK++;
			} else {

			}
		}
		K = newK;
		for (int dIndex = 0; dIndex < K; dIndex++) {
			z[dIndex] = kOldToKNew[z[dIndex]];
		}
	}
	//renumber the arr
	public void swap(int[] arr, int arg1, int arg2){
		int t = arr[arg1]; 
		arr[arg1] = arr[arg2]; 
		arr[arg2] = t; 
	}
	//renumber the arr
	public void swap(double[] arr, int arg1, int arg2){
		double t = arr[arg1]; 
		arr[arg1] = arr[arg2]; 
		arr[arg2] = t; 
	}
	//renumber the b arr
	public void swap(boolean[][] arr, int arg1, int arg2) {
		boolean[] t = arr[arg1]; 
		arr[arg1] = arr[arg2]; 
		arr[arg2] = t; 
	}
	//renumber the arr
	public  void swap(int[][] arr, int arg1, int arg2) {
		int[] t = arr[arg1]; 
		arr[arg1] = arr[arg2]; 
		arr[arg2] = t; 
	}
	//when creating a new cluster, we will ensure capacity of the arr
	public  int[] ensureCapacity(int[] arr,int i) {
		int length = arr.length;
		int[] arr2 = new int[length  + i];
		System.arraycopy(arr, 0, arr2, 0, length);
		return arr2;
	}
	//when creating a new cluster, we will ensure capacity of the arr
	public  double[] ensureCapacity(double[] arr,int i) {
		int length = arr.length;
		double[] arr2 = new double[length  + i];
		System.arraycopy(arr, 0, arr2, 0, length);
		return arr2;
	}
	//when creating a new cluster, we will ensure capacity of the arr
	public  int[][] ensureCapacity(int[][] array,int i,int j) {  
		int[][] arr = new int[array.length +i][array[0].length +j];       
		for(int c = 0; c< array.length; c++) {  
			System.arraycopy(array[c], 0, arr[c], 0, array[c].length);  
		}  
		return arr;  
	} 
	//when creating a new cluster, we will ensure capacity of the arr
	public  boolean[][] ensureCapacity(boolean[][] array,int i,int j) {  
		boolean[][] arr = new boolean[array.length +i][array[0].length +j];    
		for(int c = 0; c< array.length; c++) {  
			System.arraycopy(array[c], 0, arr[c], 0, array[c].length);  
		}  
		return arr;  
	} 
	/**
	 * update the count 
	 * 
	 * @param probs kd,nkw,nksum
	 * @return
	 */
	void updateCount(int d, int cluster, boolean flag) {
		if (flag) {  // decrease the count
			kd[cluster]--;
			for(int n = 0; n < docword[d].length; n ++) {
				nkw[cluster][docword[d][n]]--;
				nksum[cluster]--;
			}
		}else { // increase the count
			kd[cluster]++;
			for(int n = 0; n < docword[d].length; n ++) {
				nkw[cluster][docword[d][n]]++;
				nksum[cluster]++;
			}
		}
	}
	public void sampleBinaryBMatrix() {
		int GIBBS_ITER = 1;
		b_sum = new int[K];
		for (int k = 0; k != K; k++) {
			for (int v = 0; v != V; v++) {
				b[k][v] = nkw[k][v] > 0;
				b_sum[k] += b[k][v] ? 1 : 0;
			}
		}
		//
		double log_diff, ratio, p;
		for (int iter = 0; iter != GIBBS_ITER; iter++) {
			for (int k = 0; k != K; k++) {
				for (int v = 0; v != V; v++) {
					if (b[k][v] && nkw[k][v] == 0) {
						log_diff = Gamma.logGamma(b_sum[k]*beta1 + V*beta0)
								- Gamma.logGamma((b_sum[k]-1)*beta1 + V*beta0);
						log_diff -= Gamma.logGamma(nksum[k] + b_sum[k]*beta1 + V*beta0)
								- Gamma.logGamma(nksum[k] + (b_sum[k]-1)*beta1 + V*beta0);

						ratio = Math.exp(log_diff) * pi_b[k] / (1.0-pi_b[k]);
						p = ratio / (1.0 + ratio);
						if (rand.nextDouble() > p) { 
							b[k][v] = false;
							b_sum[k] --;
						}
					} else if (!b[k][v]) {
						log_diff = Gamma.logGamma((b_sum[k]+1)*beta1 + V*beta0)
								- Gamma.logGamma(b_sum[k]*beta1 + V*beta0);
						log_diff -= Gamma.logGamma(nksum[k] + (b_sum[k]+1)*beta1 + V*beta0)
								- Gamma.logGamma(nksum[k] + b_sum[k]*beta1 +  V*beta0);

						ratio = Math.exp(log_diff) * pi_b[k] / (1.0 - pi_b[k]);
						p = ratio / (1.0 + ratio);
						if (rand.nextDouble() < p) { 
							b[k][v] = true;
							b_sum[k] ++;
						}
					}
				}

				BetaDistribution betaDist = new BetaDistribution(rand, s + b_sum[k], r + V - b_sum[k]);
				pi_b[k] = betaDist.sample();
			}
		}
	}
	/**
	 * obtain the parameter Phi
	 */
	public double[][] estimatePhi() {
		double[][] phi = new double[K][V];
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < V; w++) {
				b[k][w] = nkw[k][w] > 0;
				int y = b[k][w] ? 1 : 0;
				phi[k][w] = (nkw[k][w] + y*beta1 + beta0) / (nksum[k] + V*beta0 + b_sum[k]*beta1);
			}
		}
		return phi;
	}
	/**
	 * obtain the parameter theta
	 */
	public double[] estimateTheta() {
		double[] theta = new double[K];
		for (int k = 0; k < K; k++) {
			theta[k] = (kd[k] + alpha) / (M - 1 + K * alpha);
		}
		return theta;
	}
	/**
	 * write top words with probability for each cluster
	 */
	public void writeTopWordsWithProbability(){
		StringBuilder sBuilder = new StringBuilder();
		double[][] phi = estimatePhi();
		int topicNumber = 1;
		for (double[] phi_z : phi) {
			sBuilder.append("Topic:" + topicNumber + "\n");
			for (int i = 0; i < topWordsOutputNumber; i++) {
				int max_index = FuncUtils.maxValueIndex(phi_z);
				sBuilder.append(indexToWordMap.get(max_index) + " :" + phi_z[max_index] + "\n");
				phi_z[max_index] = 0;
			}
			sBuilder.append("\n");
			topicNumber++;
		}
		try {
			FileUtil.writeFile(outputFileDirectory + "DPMM_cluster_word_" + K + ".txt", sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * write theta for each cluster
	 */
	public void writeClusterDistri(){
		double[] theta = estimateTheta();
		StringBuilder sBuilder = new StringBuilder();
		for (int k = 0; k < K; k++) {
			sBuilder.append(theta[k] + "\n");
		}
		try {
			FileUtil.writeFile(outputFileDirectory + "DPMM_theta_" + K + ".txt", sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * write cluster for each document
	 */
	public void writeCoherence(String inputfile){
		//read data
		ArrayList<String> docLines = new ArrayList<String>();
		FileUtil.readLines(inputfile, docLines,"gbk");
		int D = docLines.size();
		int[][] doc = new int[D][];
		int j = 0;
		for(String line : docLines){
			List<String> words = new ArrayList<String>();
			FileUtil.tokenizeAndLowerCase(line, words);
			doc[j] = new int[words.size()];
			for(int i = 0; i < words.size(); i++){
				String word = words.get(i);
				doc[j][i] = wordToIndexMap.get(word);
			}
			j++;
		}
		StringBuilder sBuilder = new StringBuilder();
		double[][] phi = estimatePhi();
		double average_coherence5 = EstimationUtil.average_coherence(doc,phi,5);
		double average_coherence10 = EstimationUtil.average_coherence(doc,phi,10);
		double average_coherence15 = EstimationUtil.average_coherence(doc,phi,15);
		double average_coherence20 = EstimationUtil.average_coherence(doc,phi,20);
		sBuilder.append(K + "\t Average_coherence is: \n");
		sBuilder.append("top 5:\t" + average_coherence5 + "\n");
		sBuilder.append("top 10:\t" + average_coherence10 + "\n");
		sBuilder.append("top 15:\t" + average_coherence15 + "\n");
		sBuilder.append("top 20:\t" + average_coherence20 + "\n");
		try {
			FileUtil.writeFile(outputFileDirectory + "SBDPModel_coherence" + K + ".txt", sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * write cluster for each document
	 */
	public void writeDocCluster(){
		StringBuilder sBuilder = new StringBuilder();
		for(int d = 0; d < M; d++){
			int cluster = z[d];
			sBuilder.append(cluster + "\n");
		}
		try {
			FileUtil.writeFile(outputFileDirectory + "DPMM_doc_cluster" + K + ".txt", sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void main(String[] args) {
		SBDPModel_Linux dmm = new SBDPModel_Linux("/home/qianyang/yelp/reviews_btm", "gbk", 33, 0.1,
				0.01, 2000, 10, "/home/qianyang/yelp/output/");
		dmm.MCMCSampling();
	}

}
