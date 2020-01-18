package util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GenerateBiterms {
	static Map<String, Integer> wordToIndexMap = new HashMap<String, Integer>();;  //word to index
	static List<String> indexToWordMap = new ArrayList<String>();    //index to String word 
	static int M; // number of documents in the corpus
	static int [][] docword;//word index array
	//biterm realted
	static int[][] biterms;
	static int windowSize = 100;
	static String outputfile = "reviews/reviews_btm";
	public static void main(String[] args) {
		//read data
		ArrayList<String> docLines = new ArrayList<String>();
		FileUtil.readLines("reviews/reviews_small", docLines,"gbk");
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
		biterms = generateBiterms(docword, windowSize);
		StringBuilder sBuilder = new StringBuilder();
		for(int b = 0; b < biterms.length; b++){
			sBuilder.append(indexToWordMap.get(biterms[b][0])+ " " + indexToWordMap.get(biterms[b][1]) + "\n");
		}
		try {
			FileUtil.writeFile(outputfile, sBuilder.toString(),"gbk");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * generate biterms
	 * @param documents 
	 * @param windowSize 
	 * @return biterms
	 */
	public static int[][] generateBiterms(int[][] documents, int windowSize) {
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
}
