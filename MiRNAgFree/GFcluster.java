package miRNAgFreeGit;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import libs.Util;
import libs.Write;



public class GFcluster {

	List<Integer> start;
	List<String> seqs;
	List<Integer> rc;
	String prominentSeq;// the sequence of the read with highes RC
	int prominentRC; // the RC of the prominent read
	public GFcluster(String initSeq, int prominentRC){
		this.prominentRC = prominentRC;
		this.prominentSeq = initSeq;
		this.start = new ArrayList<Integer>();
		this.seqs = new ArrayList<String>();
		this.rc = new ArrayList<Integer>();
		
	}
	
	
	public void add(int pos, String seq, int rc){
		this.start.add(pos);
		this.seqs.add(new String(seq));
		this.rc.add(rc);
	}
	
	/**
	 * gives back '-1' if cluster only has one read
	 * if return value = 0 --> no other reads start at position 0
	 * @return
	 */
	public double get5pFluct(){
		double c =  0;
		double total = 0;
		for(int i = 0; i < this.start.size() ; i++){
			total += this.rc.get(i);
			if (this.start.get(i) == 0)
				c+= this.rc.get(i);
		}
		if(total == 0)
			return -1;
		
		return c/total;
	}
	
	// calculates the total RC (read count) of the cluster
	public int getTotalRC(){
		double c = 0d;
		for(double r : this.rc)
			c+=r;
		
		c+= this.prominentRC;
		return (int)c;
	}
	
	public double getDominant2All(){
		return (double)this.prominentRC / (double)this.getTotalRC();
	}
	
	// the total number of reads in the cluster
	public int getReadNumber(){
		return this.start.size() + 1;
	}
	


	public void representCluster(String file,int upStreamOverhang){

		int i = 0;
		String fString = Util.getCharString(upStreamOverhang, ' ');

		Write.writeString(file, fString+this.prominentSeq+"\t"+this.prominentRC, false);
		for(String seq : this.seqs){

			fString = Util.getCharString(upStreamOverhang+this.start.get(i), ' ');
			Write.writeString(file, fString+seq+"\t"+this.rc.get(i), true);
			i++;
		}


	}
	
	
	public static void sort(List<GFcluster> list){
		
		
		
		Comparator<GFcluster> comparator = new Comparator<GFcluster>(){
			public int compare(GFcluster o1, GFcluster o2) {
				if (o1.getTotalRC() < o2.getTotalRC())
					return 1;
				if (o1.getTotalRC() > o2.getTotalRC())
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);
		
		
	}
	
}
