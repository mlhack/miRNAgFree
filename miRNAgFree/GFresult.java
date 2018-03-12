package miRNAgFreeGit;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;


public class GFresult {

	String knownPatternMicroRNAs;
	RNAcofoldData rnaCofoldData;
	String processing; // should be strict or lax --> depending on the Dicer/Droha pattern
	String homologName;
	String homologSequence;
	int pos5p;
	int pos3p;
	double probablility;
	public GFresult(RNAcofoldData data, String microRNAs, double probability){
		this.rnaCofoldData = data;
		this.knownPatternMicroRNAs = microRNAs;
		this.probablility = probability;
	}
	
	public GFresult(){

	}
	public int getExpression(){
		return this.rnaCofoldData.clust5p.getTotalRC();
	}
	
	/**
	 * Sort a list of predicted microRNAs (results objects) by it expression value
	 * @param list
	 */
	public static void sortBiggerToSmaller(List<GFresult> list){
		
		
		Comparator<GFresult> comparator = new Comparator<GFresult>(){
			public int compare(GFresult o1, GFresult o2) {
				if (o1.getExpression() < o2.getExpression())
					return 1;
				if (o1.getExpression() > o2.getExpression())
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);

	}
}
