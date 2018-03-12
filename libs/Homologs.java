package libs;

import java.util.Hashtable;
import java.util.Map;

public class Homologs {

	Map<String,String> seedFamMap = null;
	int seedStart = 2;
	int seedEnd = 8;
	public Homologs(String file, int seedStart, int seedEnd){
		
		this.seedStart = seedStart;
		this.seedEnd = seedEnd;
		seedFamMap = Read.getFastaMapInverseSeed(file, seedStart, seedEnd);
	}
	
	
	/**
	 * get the name of the miRNA family 
	 * @param seedSequence
	 * @return
	 */
	public String getName(String seedSequence){
		
		Map<String,int[]> counter = new Hashtable<String,int[]>();
		if(seedFamMap.containsKey(seedSequence)){
//			System.out.println("found: "+seedFamMap.get(seed));
			String[] f = seedFamMap.get(seedSequence).split("=");
			for(String ff : f){
				String famName = ff.replace("-5p", "").replace("-3p", "");
				famName = famName.substring(4, famName.length());
				if(counter.containsKey(famName)){
					counter.get(famName)[0]++;
				}	
				else{
					int[] t = new int[1];
					t[0]++;
					counter.put(famName, t);
				}
			}
			int maxC = 0;
			String back = null;
			for(String key : counter.keySet()){
				if(counter.get(key)[0] > maxC){
					maxC = counter.get(key)[0];
					back = new String(key);
				}
			}
			return back;
		}
		return null;
	}
}
