package libs;

import java.util.Map;

public class Align {

	/**
	 * 
	 * @param sequence
	 * @param adapter
	 * @param startPos the first position in the read where the adapter should be detected 0-based!!!!
	 * @param MM
	 * @param minLen
	 * @return returns the first base position of the adapter in 0-based coordinates; if not found --> returns -1
	 */
	public static int alignAdapter(String sequence, String adapter, int startPos, int MM, int minLen){
		
		int seqLen = sequence.length();
		for(int i = startPos ; i < seqLen; i++){
			
			if(i + minLen <= seqLen){
				
				int c=0;
				int m=0;
				for(int j = i; j < i+minLen; j++){
					
					if(sequence.charAt(j) != adapter.charAt(c)){
						m++;
					}
					if(m > MM){
						break;
					}
					c++;
				}
				if(m <= MM){
					return i;
				}
			}
			else{
				return -1;
			}
			
		}
		return -1;
	}

	
	/**
	 * The function searches if 'sequence' does map with maxMM mismatches to any of the sequences that are stored as the keys in map
	 * @param sequence
	 * @param map
	 * @param maxMM
	 * @return microRNA name or null
	 */
	static public String mapSeqToMap(String sequence, Map<String,String> map, int maxMM){
		
		int pos = -1;
		for(int i = 0; i <= maxMM; i++){
			for(String seq : map.keySet()){
				
				if(sequence.length() >= seq.length())
					pos = alignAdapter(sequence,seq, 0, i, seq.length());
				else{
					pos = alignAdapter(sequence,seq.substring(0, sequence.length()), 0, i, sequence.length());	
				}
				if(pos >= 0 && pos <= 1){
					return map.get(seq);
				}
			}
		}
		return null;
	}
	
	/**
	 * give back the best hit with less than 'maxMM': an array with [0] --> position; [1] --> mismatches
	 * @param query
	 * @param reference
	 * @param maxMM
	 * @return
	 */
	public static int[] alignSequence(String query, String reference, int maxMM){
		
		int stop = reference.length() - query.length();
		int qLen = query.length();
		int pos = -1;
		int mmMin = maxMM+1;
		for(int i = 0 ; i < stop ; i++){

			int m=0;
			for(int j = 0; j < qLen; j++){

				if(reference.charAt(i+j) != query.charAt(j)){
					m++;
				}
				if(m > maxMM){
					break;
				}

			}
			
			if(m < mmMin){
				pos = i;
				mmMin = m;
			}

		}

		if(pos >= 0){
			int[] back = new int[2];
			back[0] = pos;
			back[1] = mmMin;
			return back;
		}
		else
			return null;
	}
}
