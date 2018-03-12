package libs;

public class AlignData {

	public int ur; // number of unique reads mapped
	public int rc; // total number of reads mapped
	public double rcAdj; //adjusted read count
	public int urS; // unique reads mapped in sense orientation
	public int rcS; // total read count in sense orientation
	public double rcSAdj; // adjusted read count in sense orientation
	public int urAS; // unique reads mapped in anti-sense orientation
	public int rcAS; // total read count in anti-sense orientation
	public double rcASAdj; // adjusted read count in anti-sense orientation
	
	public String libraryPathBowtie;
	public String libraryPathFasta;
	public String library; // the name of the library
	public String parsedFile; // the parsed file on which the data is based
	public String senseOut; // the MA output file (sense)
	public String aSenseOut; // the SA output file (antisense)
	
	public AlignData(String parsedFile, int ur, int rc, double rcAdj, int urS, int rcS, double rcSAdj, int urAS, int rcAS, double rcASAdj){

		this.ur = ur;
		this.rc = rc;
		this.rcAdj = rcAdj;
		
		this.urS = urS;
		this.rcS = rcS;
		this.rcSAdj = rcSAdj;
		
		this.urAS = urAS;
		this.rcAS = rcAS;
		this.rcASAdj = rcASAdj;

		this.parsedFile = parsedFile;
	}
	
	
}
