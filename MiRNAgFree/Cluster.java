package miRNAgFreeGit;

import java.util.ArrayList;

import java.util.List;

import java.util.Set;

import libs.GVars;
import libs.IO;
import libs.Sort;
import libs.Util;


public class Cluster {

	

	public static List<GFcluster> cluster(String reads){
		
		/////////////////////////////////////////////////////////
		//////  prepare the reads for clustering
		/////   i) get the reads sorted from high to low copy number
		/////   ii) remove the last base (might be a NTA)
		
		List<Sort> list = Util.getSRNAbenchFormatSortList(reads);
		List<String> removedBasesList = GfreeHelper.getLengthMinusOne(list,1);
		if(list.size() != removedBasesList.size()){
			IO.warning("Error after reads preparation. Will quite now");
			IO.log(GVars.logFile, 4, "Error after reads preparation. Will quite now", true);
			System.exit(1);
		}
		IO.writeToCommandLineL1("Start analysis with "+list.size()+" unique reads");
		IO.log(GVars.logFile, 2, "Start analysis with "+list.size()+" unique reads", true);	
		
		//////////////////////////////////////////////////////////
		//////////////
		////////////////////////////////////////////////////////////
		
		List<GFcluster> back = new ArrayList<GFcluster>();
		for(int i = 0; i < list.size(); i++){

			GFcluster clust = new GFcluster(list.get(i).line,(int)list.get(i).sortValue);
			getClusterNotEqual(list, removedBasesList.get(i), clust,removedBasesList);
			back.add(clust);
			list.remove(i);
			removedBasesList.remove(i);
			i--;
		}
		return back;
	}
	
	/**
	 * the same function as getCluster however this functions does not check for equal --> speed up
	 * @param list
	 * @param seq
	 * @param cluster
	 */
	public static void getClusterNotEqual(List<Sort>list, String seq, GFcluster cluster, List<String> removedBasesList){
		
		for(int i = 1; i < list.size() ; i++){
			int pos = alignTwoSequencesMinLen(seq, removedBasesList.get(i), Vars.mm, Vars.overhang5p, Vars.overhang3p);
			

//			System.out.println(pos+"  "+ seq5p+"\n"+list.get(i).line);
			if(pos != -99){
//				System.out.println(pos+" "+list.get(i).sortValue+"  \n"+ seq+"\n"+list.get(i).line);
				cluster.add(pos, list.get(i).line, (int)list.get(i).sortValue);
				list.remove(i);
				removedBasesList.remove(i);
				i--;			
			}
		}
	}
	
	
	/**
	 * Align two sequences with 3 bp fluctuation at the 5p 
	 * @param seq1
	 * @param seq2
	 * @param MM
	 * @param seq2Overhang - number of bases that the seq2 can overhang at 5p (number of bases that do not align to seq1 at 5p)
	 * @return
	 */
	public static int alignTwoSequencesMinLen(String seq1, String seq2, int MM, int seq2Overhang, int seq1Overhang){
		
		int seqLen1 = seq1.length();
		int seqLen2 = seq2.length();
		// check 5p overhang
		for(int i = 0 ; i <= seq2Overhang; i++){

			int m=0; // mismatch count
			for(int j = 0; j < seqLen1; j++){
				
				if(j >= seqLen1 || j+i >= seqLen2){
						break;
				}
				
				if(seq1.charAt(j) != seq2.charAt(j+i)){
					m++;
				}
				if(m > MM){
					break;
				}

			}
			if(m <= MM ){
				return -i;
			}
		}

		for(int i = 0 ; i <= seq1Overhang; i++){
			int m=0;
			for(int j = 0; j < seqLen2; j++){

				if(j+i >= seqLen1 || j >= seqLen2){
					break;
				}

				if(seq1.charAt(j+i) != seq2.charAt(j)){
					m++;
				}
				if(m > MM){
					break;
				}
			}
			if(m <= MM ){
				return i;
			}
		}

		return -99;
	}
	
	
	/**
	 * go from the last to the first entry and remove all Sort objects from the list that have been assigned
	 * @param list
	 * @param detected
	 */
	public static void removeDetected(List<Sort>list, Set<Sort> detected){
		int length = list.size();
		for(int i = length -1; i >= 0; i--){
			if(detected.contains(list.get(i))){
				list.remove(i);
			}
		}
	}
	
	
}
