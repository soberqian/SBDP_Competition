package model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import evaluation.EstimationUtil;
import util.FileUtil;
import util.FuncUtils;

/**
 * 
 * Collapsed Gibbs sampling in the generative model of Latent Dirichlet Allocation
 * 
 * Reference:
 * Griffiths T. Gibbs sampling in the generative model of latent dirichlet allocation[J]. 2002.
 * Heinrich G. Parameter estimation for text analysis[R]. Technical report, 2005.
 * 
 * @author: Yang Qian
 */

public class GibbsSamplingLDA
{
	public double alpha; // Hyper-parameter alpha
	public double beta; // Hyper-parameter beta
	public int K; // number of topics
	public int iterations; // number of Gibbs sampling iterations
	public int topWords; // number of most probable words for each topic
	public Map<String, Integer> wordToIndexMap = new HashMap<String, Integer>();;  //word to index
	public List<String> indexToWordMap = new ArrayList<String>();    //index to String word 
	public int [][] docword;//word index array
	public int M; // number of documents in the corpus
	public int V; // number of words in the corpus
	public int[][] ndk; // document-topic count
	public int[] ndsum; //document-topic sum
	public int[][] nkw; //topic-word count
	public int[] nksum; //topic-word sum (total number of words assigned to a topic)
	public int[][] z;
	//output
	public int topWordsOutputNumber;
	public String outputFileDirectory; 
	public GibbsSamplingLDA(String inputFile, String inputFileCode, int topicNumber,
			double inputAlpha, double inputBeta, int inputIterations, int inTopWords,
			String outputFileDir){
		//read data
		ArrayList<String> docLines = new ArrayList<String>();
		FileUtil.readLines(inputFile, docLines,inputFileCode);
		M = docLines.size();
		docword = new int[M][];
		int j = 0;
		for(String line : docLines){
			List<String> words = new ArrayList<String>();
			FileUtil.tokenizeAndLowerCase(line, words);
			docword[j] = new int[words.size()];
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
			}
			j++;

		}
		V = indexToWordMap.size();
		alpha = inputAlpha;
		beta = inputBeta;
		K = topicNumber;
		iterations = inputIterations;
		topWordsOutputNumber = inTopWords;
		outputFileDirectory = outputFileDir;
		initialize();
	}
	/**
	 * Randomly initialize topic assignments
	 */
	public void initialize(){
		ndk = new int[M][K];
		ndsum = new int[M];
		nkw = new int[K][V];
		nksum = new int[K];
		z = new int[M][];
		for (int d = 0; d < M; d++) {
			int Nd = docword[d].length;  // the number of words in a document
			z[d] = new int[Nd];
			for (int n = 0; n < Nd; n++) {
				int topic = (int) (Math.random() * K);
				z[d][n] = topic;
				updateCount(d, topic, docword[d][n], +1);
			}
		}
	}
	public void MCMCSampling(){
		for (int iter = 1; iter <= iterations; iter++) {
			System.out.println("iteration : " + iter);
			gibbsOneIteration();
		}
		writeCoherence();
	}
	public void gibbsOneIteration() {
		for (int d = 0; d < M; d++) {
			for (int n = 0; n < z[d].length; n++) {
				int topic = z[d][n]; // get the old topic
				updateCount(d, topic, docword[d][n], -1); // update the count --1
				double[] p = new double[K];
				for (int k = 0; k < K; k++) {
					p[k] = (ndk[d][k] + alpha) / (ndsum[d] + K * alpha) * (nkw[k][docword[d][n]] + beta)
							/ (nksum[k] + V * beta);
				}
				topic = FuncUtils.rouletteGambling(p); //roulette gambling for updating the topic of a word
				z[d][n] = topic;
				updateCount(d, topic, docword[d][n], +1);  // update the count ++1

			}
		}
	}
	/**
	 * update the count 
	 * 
	 * @param probs ndk,ndsum,nkw,nksum
	 * @return
	 */
	void updateCount(int d, int topic, int word, int flag) {
		ndk[d][topic] += flag;
		ndsum[d] += flag;
		nkw[topic][word] += flag;
		nksum[topic] += flag;
	}
	/**
	 * obtain the parameter Phi
	 */
	public double[][] estimatePhi() {
		double[][] phi = new double[K][V];
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < V; w++) {
				phi[k][w] = (nkw[k][w] + beta) / (nksum[k] + V * beta);
			}
		}
		return phi;
	}
	/**
	 * write cluster for each document
	 */
	public void writeCoherence(){
		StringBuilder sBuilder = new StringBuilder();
		double[][] phi = estimatePhi();
		double average_coherence5 = EstimationUtil.average_coherence(docword,phi,5);
		double average_coherence10 = EstimationUtil.average_coherence(docword,phi,10);
		double average_coherence15 = EstimationUtil.average_coherence(docword,phi,15);
		double average_coherence20 = EstimationUtil.average_coherence(docword,phi,20);
		sBuilder.append(K + "\t Average_coherence is: \n");
		sBuilder.append("top 5:\t" + average_coherence5 + "\n");
		sBuilder.append("top 10:\t" + average_coherence10 + "\n");
		sBuilder.append("top 15:\t" + average_coherence15 + "\n");
		sBuilder.append("top 20:\t" + average_coherence20 + "\n");
		try {
			FileUtil.writeFile(outputFileDirectory + "LDA_coherence" + K + ".txt", sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void main(String args[]) throws Exception{
		GibbsSamplingLDA lda = new GibbsSamplingLDA("/home/qianyang/yelp/reviews_small", "gbk", 24, 0.1,
				0.01, 2000, 50, "/home/qianyang/yelp/output/");
		lda.MCMCSampling();
	}
}
