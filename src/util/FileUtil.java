package util;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
public class FileUtil {
	//read a file to list
	public static void readLines(String file, ArrayList<String> lines, String code) {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader( new InputStreamReader( new FileInputStream( new File(file)),code));
			String line = null;
			while ((line = reader.readLine()) != null) {
				lines.add(line);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	// write list to a file
	public static void writeLines(String file, ArrayList<?> counts, String code) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter( new OutputStreamWriter( new FileOutputStream( new File(file)),code));
			for (int i = 0; i < counts.size(); i++) {
				writer.write(counts.get(i) + "\n");
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	/**
	 * 
	 * @param file
	 * @param content
	 * @param code
	 * @throws IOException
	 */
	public static void writeFile(String file, String content,String code) throws IOException {

		File fileOutput = new File(file);
		OutputStream out = new FileOutputStream(fileOutput, false);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out, code));
		bw.write(content);
		bw.close();
		out.close();
	}
	/**
	 * split the document
	 * @param line
	 * @param tokens
	 * */
	public static void tokenizeAndLowerCase(String line, List<String> tokens) {
		StringTokenizer strTok = new StringTokenizer(line);
		while (strTok.hasMoreTokens()) {
			String token = strTok.nextToken();
			tokens.add(token.toLowerCase().trim());
		}
	}
}
