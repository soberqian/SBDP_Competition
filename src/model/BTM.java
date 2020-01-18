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
 * Gibbs sampling for BTM
 * 
 * Reference:
 * Cheng X, Yan X, Lan Y, et al. Btm: Topic modeling over short texts[J]. IEEE Transactions on Knowledge and Data Engineering, 2014, 26(12): 2928-2941.
 * Yan X, Guo J, Lan Y, et al. A biterm topic model for short texts[C]//Proceedings of the 22nd international conference on World Wide Web. ACM, 2013: 1445-1456.
 * 
 * @author: Yang Qian
 */

public class BTM
{
	public double alpha; // Hyper-parameter alpha
	public double beta; // Hyper-parameter beta
	public int K; // number of topics
	public int iterations; // number of iterations
	public Map<String, Integer> wordToIndexMap = new HashMap<String, Integer>();;  //word to index
	public List<String> indexToWordMap = new ArrayList<String>();    //index to String word 
	public int M; // number of documents in the corpus
	public int V; // number of words in the corpus
	public int [][] docword;//word index array
	//biterm realted
	public int[][] biterms;
	public int windowSize;
	public int[] z;
	public int[][] nkw; 
	public int[] nkw_sum; 
	public int[] nk;
	//output
	public int topWordsOutputNumber;
	public String outputFileDirectory; 
	public BTM(String inputFile, String inputFileCode, int topicNumber,
			double inputAlpha, double inputBeta, int inputIterations, int inTopWords, int windowS,
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
		windowSize = windowS;
		iterations = inputIterations;
		topWordsOutputNumber = inTopWords;
		outputFileDirectory = outputFileDir;
		//generate biterms
		biterms = generateBiterms(docword, windowSize);
		//initialize
		initialize();
	}
	/**
	 * Randomly assign the topic for each biterm
	 */
	public void initialize(){
		//Biterm size
		int NB = biterms.length;
		//biterm realted
		z = new int[NB];
		nkw = new int[K][V];
		nkw_sum = new int[K];
		nk = new int[K];
		for (int b = 0; b < NB; ++b) {
			int topic = (int) (Math.random() * K);
			z[b] = topic;
			nkw[topic][biterms[b][0]]++;
			nkw[topic][biterms[b][1]]++;
			nk[topic]++;
			nkw_sum[topic] += 2;
		}
	}
	public void MCMCSampling(){
		for (int iter = 1; iter <= iterations; iter++) {
			System.out.println("iteration : " + iter);
			gibbsOneIteration();
		}
		// output the result
		writeCoherence();

	}
	public void gibbsOneIteration() {
		for (int i = 0; i < biterms.length; i++) {
			int topic = z[i];
			updateCount(i, topic, 0);
			double[] p = new double[K];
			for (int k = 0; k < K; ++k) {
				p[k] = (nk[k] + alpha) * ((nkw[k][biterms[i][0]] + beta) / (nkw_sum[k] + V * beta))
						* ((nkw[k][biterms[i][1]] + beta) / (nkw_sum[k] + V * beta + 1));

			}
			topic = FuncUtils.rouletteGambling(p); //roulette gambling for updating the topic of a word
			z[i] = topic;
			updateCount(i, topic, 1);
		}
	}
	/**
	 * update the count nkw, nk and nkw_sum
	 * 
	 * @param biterm
	 * @param topic
	 * @param flag
	 * @return null
	 */
	void updateCount(int biterm, int topic, int flag) {
		if (flag == 0) {
			nkw[topic][biterms[biterm][0]]--;
			nkw[topic][biterms[biterm][1]]--;
			nk[topic]--;
			nkw_sum[topic] -= 2;
		}else {
			nkw[topic][biterms[biterm][0]]++;
			nkw[topic][biterms[biterm][1]]++;
			nk[topic]++;
			nkw_sum[topic] += 2;
		}

	}
	/**
	 * obtain the parameter phi
	 */
	public double[][] estimatePhi() {
		double[][] phi = new double[K][V];
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < V; w++) {
				phi[k][w] = (nkw[k][w] + beta) / (nkw_sum[k] + V * beta);
			}
		}
		return phi;
	}
	/**
	 * generate biterms
	 * @param documents 
	 * @param windowSize 
	 * @return biterms
	 */
	public int[][] generateBiterms(int[][] documents, int windowSize) {
		List<int[]> list = new ArrayList<>();
		for (int d = 0; d < documents.length; ++d) {
			for (int i = 0; i < documents[d].length - 1; ++i) {
				for (int j = i + 1; j < Math.min(i + windowSize, documents[d].length); ++j) {
					list.add(new int[]{documents[d][i], documents[d][j]});
				}
			}
		}
		int[][] biterms = new int[list.size()][2];
		list.toArray(biterms);
		return biterms;
	}
	/**
	 * generate biterms for a document
	 * @param documents 
	 * @param windowSize 
	 * @return biterms
	 */
	public int[][] generateBitermsForOneDoc(int[] document, int windowSize) {
		List<int[]> list = new ArrayList<>();
		for (int i = 0; i < document.length - 1; ++i) {
			for (int j = i + 1; j < Math.min(i + windowSize, document.length); ++j) {
				list.add(new int[]{document[i], document[j]});
			}
		}
		int[][] biterms = new int[list.size()][2];
		list.toArray(biterms);
		return biterms;
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
			FileUtil.writeFile(outputFileDirectory + "BTM_coherence" + K + ".txt", sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void main(String args[]) throws Exception{
		BTM btm = new BTM("/home/qianyang/yelp/reviews_small", "gbk", 35, 0.1,
				0.01, 2000, 50, 50, "/home/qianyang/yelp/output/");
		btm.MCMCSampling();
	}
}
