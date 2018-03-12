package miRNAgFreeGit;

public class AlignData {

	int position;
	String name;
	String sequence;
	int alignLength;
	int mm;
	/**
	 * Holds the alignment information of a read against a reference sequence
	 * @param name of the sequence
	 * @param sequence seq. of reference sequence (i.e. putative homologous)
	 * @param mm number of mismatches
	 * @param pos position of the alignment
	 * @param alignLength  length of the alignment
	 */
	public AlignData(String name, String sequence, int mm, int pos, int alignLength){
		
		this.name = name;
		this.sequence = sequence;
		this.mm = mm;
		this.position = pos;
		this.alignLength = alignLength;
		
	}
}
