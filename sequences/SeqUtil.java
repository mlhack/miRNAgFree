package sequences;

import java.util.Hashtable;
import java.util.Map;


/**
 * Functions for (short) sequence manipulation
 * @author Michael
 *
 */
public class SeqUtil {

	private static Map<Character, Character> revbas;

	static{
		revbas = new Hashtable<Character,Character>();
		revbas.put('A', 'T');
		revbas.put('C', 'G');
		revbas.put('G', 'C');
		revbas.put('T', 'A');
		revbas.put('-', '-');
		revbas.put('N', 'N');
		revbas.put('a', 't');
		revbas.put('c', 'g');
		revbas.put('g', 'c');
		revbas.put('t', 'a');
		revbas.put('n', 'n');
		revbas.put('R', 'Y');
		revbas.put('r', 'y');
		revbas.put('Y', 'R');
		revbas.put('y', 'r');
		
		revbas.put('M', 'K');
		revbas.put('K', 'M');
		revbas.put('S', 'S');
		revbas.put('W', 'W');
		
		revbas.put('H', 'D');
		revbas.put('D', 'H');
		
		revbas.put('B', 'V');
		revbas.put('V', 'B');
		
	}

	/**
	 * The function gives back the reverse complementary of the input sequence
	 * @param sequence
	 * @return
	 */
	public static String getReverseComplementarySequence(String sequence){

		StringBuilder newSeq = new StringBuilder(sequence.length());
		char[] seqArr = sequence.toCharArray();
		int start = seqArr.length;
		for(int i = start-1; i >= 0;i--){
			newSeq.append(revbas.get(seqArr[i]));
		}
		return newSeq.toString();
	}
	
	/**
	 * The function gives back the reverse (not complementary!!!!) of the input sequence
	 * @param sequence
	 * @return
	 */
	public static String getReverseSequence(String sequence){

		StringBuilder newSeq = new StringBuilder(sequence.length());
		char[] seqArr = sequence.toCharArray();
		int start = seqArr.length;
		for(int i = start-1; i >= 0;i--){
			newSeq.append(seqArr[i]);
		}
		return newSeq.toString();
	}
	
	/**
	 * gets the complementary but not reverse!!! of the input seq
	 * @return
	 */
	public static String getComplementarySequence(String sequence){

		StringBuilder newSeq = new StringBuilder(sequence.length());
		char[] seqArr = sequence.toCharArray();
		int start = seqArr.length;
		for(int i = 0; i < start;i++){
			newSeq.append(revbas.get(seqArr[i]));
		}
		return newSeq.toString();

	}

	


	/**
	 * gives back the ratio of the most frequent nucleotide
	 * @param read
	 * @return
	 */
	public static double ratioOfMostFrequentBase(String read){

		int[] freq = new int[4];
		int total = 0;
		for(int i = 0; i < read.length();i++){
			total++;
			if(read.substring(i, i+1).equals("A")){
				freq[0]++;
			}
			else if(read.substring(i, i+1).equals("C")){
				freq[1]++;
			}
			else if(read.substring(i, i+1).equals("G")){
				freq[2]++;
			}
			else if(read.substring(i, i+1).equals("T")){
				freq[3]++;
			}
		}

		int most = 0;
		for(int i = 0; i < freq.length; i++){
			if(freq[i] > most)
				most = freq[i];
		}
		double ratio = most/(double)total;
		return ratio;


	}
	
	
	/**
	 * Counts the number of occurencies of a pattern in a given sequence 
	 * @param sequence The sequence String 
	 * @param pattern The pattern String
	 * @return the number of occurencies int
	 */
	static public int count(String sequence, String pattern){

		int count = 0;
		int stop = sequence.length() - pattern.length() + 1;
		for(int i = 0 ; i < stop; i++){

			String sub = sequence.substring(i, i+pattern.length());
			if(sub.equalsIgnoreCase(pattern)){
				count++;
			}
		}
		return count;
	}
	

}
