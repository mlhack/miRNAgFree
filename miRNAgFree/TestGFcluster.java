package miRNAgFreeGit;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import libs.GVars;
import libs.IO;
import libs.Write;

/**
 * Takes a list of GFclusters and one GFcluster object and checks if the main object is the guide strand and one of the objects in the list
 * can be the passenger strand
 * @author michael
 *
 */
public class TestGFcluster implements Runnable{

	List<GFcluster> clusterList;
	int index;
	String tmp;
	GFresult result;
	
	public TestGFcluster(List<GFcluster> clusterList, int index, String tmp){
		this.clusterList = clusterList;
		this.index = index;
		this.tmp = tmp;
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		
		String homoName = null;
		// Vars.homologMapMature is initialized if microRNA= was given at the input
		if(Vars.homologMapMature != null){
			// if the seed sequence of the dominant read is contained in homologMapMature
			// --> assign the corresponding name
			homoName = Vars.homologs.getName(clusterList.get(this.index).prominentSeq.substring(Vars.seedStart-1, Vars.seedEnd));
					
			// if a name was assigned
			if(homoName != null){
				// test if a pre-microRNA sequence exist to which this read cluster maps & a passenger read cluster
				result = testHomolog(clusterList, index,Vars.guideMM, Vars.passMM);
				if(result != null){
					result.processing = "strict homologous";
//					IO.writeToCommandLineL2("index: "+index+" homologous: "+result.homologName+":"+result.pos5p+"-"+result.pos3p+"\n");

				}
			}
		}
		// if no homolog was found 
		if(result == null){
			// write the input file for RNAcofold
			writeTempSequencesOne();
			// launch RNAcofold and parse it
			List<RNAcofoldData> list = getHybrids(tmp);
			
			
			if(list.size() > 0){
				// get the passenger read cluster (that forms a duplex with this read cluster)
				// that has the best energy ratio & bindings above threshold & Dicer/Drosha criterion
				result = getResultML(list, clusterList, homoName,tmp);
				if(result != null){
				result.processing = "duplex";
				}
			}
		}

		// check for simple homologs (only guide sequence is in data set --> no passenger was found)
		if(result == null && homoName != null){
//			IO.warning("Found simple homologous candidate (only the guide strand is found in homologous");
			IO.log(GVars.logFile, 1, "Found lax homologous", true);
			GFresult resH = new GFresult();
			resH.homologName = homoName;
			resH.processing = "lax homologous";
			RNAcofoldData rna = new RNAcofoldData();
			rna.clust5p = clusterList.get(index);
			rna.name = index+"#-1#-1#1";
		}
		
		new File(tmp).delete();
	}
	
	
	/**
	 * gets a GFresults object based on machine learning
	 * @param list
	 * @param clusterList
	 * @param homoName
	 * @param i
	 * @return
	 */
	public  GFresult getResultML(List<RNAcofoldData> list, List<GFcluster> clusterList, String homoName, String arffFile){
		

		if(list.size() == 0)
			return null;
		
		
		RNAcofoldData.sort(list);


		if(Vars.perfectHybrid){
			if(list.get(0).bindings >= Vars.bindings){
				setClusterData(list.get(0), clusterList);

				GFresult res1 = new GFresult(list.get(0),null, 1);
				res1.homologName = homoName;
				return res1;
			}
			else{
				return null;
			}
			
		}
		else{
			for(int i = 0; i < list.size(); i++){

				if(list.get(i).bindings >= Vars.bindings && list.get(i).hasInnerNonBindings()){
					setClusterData(list.get(i), clusterList);

					GFresult res1 = new GFresult(list.get(i),null, 1);
					res1.homologName = homoName;
					return res1;
				}

			}
		}

		return null;
//
	}



	
	/**
	 * leaves only the RNAcofoldData objects in the list if the drosha-dicer fluctuation is equal or below the threshold
	 * @param cofoldList
	 * @param maxFluc
	 */
	public static void getBestFluc(List<RNAcofoldData> cofoldList, int maxFluc){
		
		for(int i = 0; i < cofoldList.size(); i++){
			if(cofoldList.get(i).droshaDicerFluc > maxFluc){
				cofoldList.remove(i);
				i--;
			}
		}
	}
	
	


	

	

	/**
	 * this function runs RNAcofold on the Helper.tmp file and gives back the "best" hybrid as a RNAcofoldData object
	 * @return
	 */
	private static List<RNAcofoldData> getHybrids(String file){

		List<String> duplex = RNAcofoldData.exec(GVars.RNAcofold, GVars.RNAcofoldParameter, file);
		List<RNAcofoldData> cofoldList  = null;
		if(Vars.procMode.equals("strict"))
			 cofoldList = RNAcofoldData.getCofoldDataStrict(duplex);
		else if(Vars.procMode.equals("lax"))
			 cofoldList = RNAcofoldData.getCofoldDataLax(duplex);
		else {  // accept all hybrids with more than 8 bindings
			cofoldList = RNAcofoldData.getCofoldDataMinBinding(duplex, 8);
		}
		return cofoldList;
		

	}

	


	public void writeTempSequencesOne(){
		/////////////////////////////////////////////////////////////////////
		//make the first temporary file only with the predominent sequence
		makeTempFileClusterOne(clusterList, index, tmp,false);
	}


	

	/**
	 * 	 * make a temprary file 
	 * @param sequence
	 * @param list
	 * @param file
	 * @param bannIndex --> the current index in the clusterList object
	 */
	public static void makeTempFileClusterOne(List<GFcluster> list, int mainIndex,  String outFile, boolean append){

		try {
			BufferedWriter writer = new BufferedWriter ( new FileWriter(outFile, append));

			for(int i = mainIndex+1; i < list.size(); i ++){

				
				// the id has the following order: [0]-> index of higher expressed cluster (initially assumed to be 5p)
				// [1] -> the index in the cluster (-1 = dominant read, other are the index in the seq list)
				// [2] --> index of 3p cluster in clusterList 
				// [3] --> index of read in 3p Clusterobject (-1 = dominant read, other are the index in the seq list)
				writer.write(">"+mainIndex+"#"+-1+"#"+i+"#"+-1+"#5p"+"\n");
				writer.write(list.get(mainIndex).prominentSeq+"&"+list.get(i).prominentSeq+"\n");
								
				for(int j = 0; j < list.get(i).seqs.size(); j++){
					
					//use only a limited number of reads
					if(j >= Vars.nrNonDomStar){
						break;
					}
					writer.write(">"+mainIndex+"#"+-1+"#"+i+"#"+j+"#5p"+"\n");
					writer.write(list.get(mainIndex).prominentSeq+"&"+list.get(i).seqs.get(j)+"\n");	
					
				}
				
				///////////////////////////////////////////////
				///// 
				////   
				for(int ii = 0; ii < Vars.nrNonDomGuide; ii++){
					
					writer.write(">"+mainIndex+"#"+ii+"#"+i+"#"+-1+"#5p\n");
					writer.write(list.get(mainIndex).seqs.get(ii)+"&"+list.get(i).prominentSeq+"\n");
					
					
					for(int j = 0; j < list.get(i).seqs.size(); j++){
						
						//use only a limited number of reads
						if(j >= Vars.nrNonDomStar){
							break;
						}
						writer.write(">"+mainIndex+"#"+ii+"#"+i+"#"+j+"#5p\n");
						writer.write(list.get(mainIndex).seqs.get(ii)+"&"+list.get(i).seqs.get(j)+"\n");			
						
					}
				}

			}
			

			
			
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	



	/**
	 * This function takes the first GFcluster object (index 0 == guide strand) and searches for a homologous microRNA
	 * If it findes one, it looks for alignments of the dominant reads of all other clusters
	 * @param list
	 * @return
	 */
	public  GFresult testHomolog( List<GFcluster> clusterList,  int mainIndex, int mainMM, int starMM){


		List<AlignData> detectedHomologous = new ArrayList<AlignData>();
		GFcluster mainCluster = clusterList.get(mainIndex);
		/////////////////////////////////////////////////
		////////
		/////// try to align the prominent sequence of the main cluster against the homologous sequences

		for(String homologSequence : Vars.homologMapHairpin.keySet()){

			int[] align = libs.Align.alignSequence(mainCluster.prominentSeq.substring(0, mainCluster.prominentSeq.length()-1), 
					homologSequence, mainMM); // align the 

			if(align != null){
				detectedHomologous.add(new AlignData(Vars.homologMapHairpin.get(homologSequence), homologSequence,
						align[1],align[0],mainCluster.prominentSeq.length()-1));
			}
			else {
				//				System.out.println("not detected "+homologMapHairpin.get(homologSequence));
			}
		}
		/////////////////////////////////////////////
		// get the best alignments of the putative guide strand --> delete the rest
		getBestAlignData(detectedHomologous);


		if(detectedHomologous.size() == 0)
			return null;


		//////////////////////////////////////////////////////////////////////
		///////////
		//////    align the putative passenger strand only to those selected in the first step

		List<Integer> bestStar = new ArrayList<Integer>(); // the indexes (in the clusterList) of the best putative star sequences are recorded
		List<Integer> bestPos = new ArrayList<Integer>(); // the positions of the best putative star sequences are recorded
		List<Integer> guideIndex = new ArrayList<Integer>(); // the corresponding index in the 'detectedHomologous' list
		int bestMM = 100;
		
		/// go over best candidates 
		for(int j = 0; j < detectedHomologous.size(); j++){

			for(int i = mainIndex+1; i < clusterList.size(); i++){


				int[] align = libs.Align.alignSequence(clusterList.get(i).prominentSeq, 
						detectedHomologous.get(j).sequence, starMM); // align the 

				//////////////////////////////////////////
				// less number of indexes --> remove previous and add new one
				if(align != null && align[1] < bestMM ){
				
					bestStar = new ArrayList<Integer>();
					bestPos = new ArrayList<Integer>();
					guideIndex = new ArrayList<Integer>();
					bestStar.add(i);
					bestPos.add(align[0]);
					guideIndex.add(j);
					bestMM = align[1];
					
				}
				// same number of mismatches --> add
				else if(align != null && align[1] == bestMM ){
					bestStar.add(i);
					bestPos.add(align[0]);
					guideIndex.add(j);
				}
			}
		}
		
		if(bestStar.size() > 0){
			return getBestHomologousHybrid(bestStar, bestPos,guideIndex, clusterList,mainIndex,detectedHomologous);
		}
		else{
			IO.log(GVars.logFile, 3, "Found guide but not passenger for index "+mainIndex, true);
			return null;
		}

	}
	
	/**
	 * Get best hybrid from all possible combinations (guide and passenger mappings). Criterion: i) the Drosha/Dicer fluctuation, ii) the total read count of the cluster
	 * @param bestStar
	 * @param bestPos
	 * @param guideIndex
	 * @param clusterList
	 * @param mainIndex
	 * @param detectedHomologous
	 * @return
	 */
	public  GFresult getBestHomologousHybrid(List<Integer> bestStar, List<Integer> bestPos, List<Integer> guideIndex,  
			List<GFcluster> clusterList,  int mainIndex,List<AlignData> detectedHomologous ){

		int bestI = -1; // the best index;
		int bestPassangerIndex = 10000000;
		RNAcofoldData best = null; // the best hybrid
		int bestFluc = 100;

		
		for(int i = 0; i < bestStar.size(); i++){
			
			int passengerIndex = bestStar.get(i);
			int passengerPos = bestPos.get(i);
			int gIndex = guideIndex.get(i);
			
			/// detect which one is 5p/3p
			/////////////////////////////////
			//  the guide is 3p


			if(passengerPos < detectedHomologous.get(gIndex).position){

				String header = ">"+passengerIndex+"#-1#"+mainIndex+"#-1";	
				Write.writeFasta(tmp, clusterList.get(passengerIndex).prominentSeq+"&"+clusterList.get(mainIndex).prominentSeq, header, "", false);

			}
			//////////////////////////////
			//// the guide is 5p
			else{

				String header = ">"+mainIndex+"#-1#"+passengerIndex+"#-1";	
				Write.writeFasta(tmp, clusterList.get(mainIndex).prominentSeq+"&"+clusterList.get(passengerIndex).prominentSeq, header, "", false);
			}
			
			List<String> duplex = RNAcofoldData.exec(GVars.RNAcofold, GVars.RNAcofoldParameter, tmp);
			List<RNAcofoldData> cofoldList = RNAcofoldData.getCofoldDataMinBinding(duplex, 0);
			if(cofoldList.size() > 0){
				int fluc = cofoldList.get(0).getDroshaDicerFluctuation();
				
				if(fluc < bestFluc){
					bestFluc = fluc;
					best = cofoldList.get(0);
					bestI = i;
					bestPassangerIndex = passengerIndex;
				}
				else if(fluc == bestFluc && passengerIndex < bestPassangerIndex){
					best = cofoldList.get(0);
					bestI = i;
					bestPassangerIndex = passengerIndex;
				}
			}
		}
		
		if(bestFluc >= 4)
			return null;


		if(bestI >= 0){
						
			setClusterData5p3p(best, clusterList);
			GFresult res = new GFresult(best,"---",-1);
			res.homologName = new String(detectedHomologous.get(guideIndex.get(bestI)).name);
			res.homologSequence = new String (detectedHomologous.get(guideIndex.get(bestI)).sequence);
			String[] ids = res.rnaCofoldData.name.trim().split("#");
			
			// the index of 5p is smaller than of 3p --> the guide strand is 5p
			if(Integer.parseInt(ids[0])  < Integer.parseInt(ids[2])){
				res.pos5p = detectedHomologous.get(guideIndex.get(bestI)).position;
				res.pos3p = bestPos.get(bestI);
			}
			else{
				res.pos3p = detectedHomologous.get(guideIndex.get(bestI)).position;
				res.pos5p = bestPos.get(bestI);
			}

			return res;
		}
		else{
			IO.warning("Hybrid Structure did not fit for putative homologous "+detectedHomologous.get(0).name);
			IO.log(GVars.logFile, 3, "Hybrid Structure did not fit for putative homologous "+detectedHomologous.get(0).name, true);
			return null;
		}
	}

	/**
	 * 
	 * @param list
	 */
	public static void getBestAlignData(List<AlignData> list){

		int bestMM = 100;
		for(AlignData data : list){
			if(data.mm < bestMM){
				bestMM = data.mm;
			}
		}

		for(int i = 0; i < list.size(); i++){
			if(list.get(i).mm > bestMM){
				list.remove(i);
				i--;
			}
		}
	}

	/**
	 * the originally asumed arm was correct
	 */
		
	public void setClusterData(RNAcofoldData data, List<GFcluster> clusterList){
		// the id has the following order: [0]-> index of more expressed cluster (originally asumed to be 5p)
		// [1] -> the index in the cluster (-1 = dominant read, other are the index in the seq list)
		// [2] --> index of 3p cluster in clusterList 
		// [3] --> index of read in 3p Clusterobject (-1 = dominant read, other are the index in the seq list)
		String[] ids = data.name.split("#");
		data.clust5p = clusterList.get(Integer.parseInt(ids[0]));
		data.clust3p = clusterList.get(Integer.parseInt(ids[2]));
	}
	public void setClusterData5p3p(RNAcofoldData data, List<GFcluster> clusterList){
		// the id has the following order: [0]-> index of more expressed cluster (originally asumed to be 5p)
		// [1] -> the index in the cluster (-1 = dominant read, other are the index in the seq list)
		// [2] --> index of 3p cluster in clusterList 
		// [3] --> index of read in 3p Clusterobject (-1 = dominant read, other are the index in the seq list)
		String[] ids = data.name.split("#");
		data.clust5p = clusterList.get(Integer.parseInt(ids[0]));
		data.clust3p = clusterList.get(Integer.parseInt(ids[2]));
	}
	/**
	 * the originally asumed arm was NOT correct --> reverse
	 */
	public  RNAcofoldData setClusterData3p5p(RNAcofoldData data, List<GFcluster> clusterList){
		
		// the id has the following order: [0]-> index of more expressed cluster (originally asumed to be 5p)
		// [1] -> the index in the cluster (-1 = dominant read, other are the index in the seq list)
		// [2] --> index of 3p cluster in clusterList 
		// [3] --> index of read in 3p Clusterobject (-1 = dominant read, other are the index in the seq list)
		
		String[] ids = data.name.split("#");
		String newID = ids[2]+"#"+ids[3]+"#"+ids[0]+"#"+ids[1];
		RNAcofoldData d = reverse(data, tmp);
		d.clust5p = clusterList.get(Integer.parseInt(ids[2]));
		d.clust3p = clusterList.get(Integer.parseInt(ids[0]));
		d.name  = newID;
		return d;
	}

	public static RNAcofoldData reverse(RNAcofoldData d, String file){
		
		Write.writeString(file, ">"+d.name, false);
		Write.writeString(file, d.sec3p+"&"+d.sec5p, true);
		RNAcofoldData back = getHybrid(file);
		return back;
	}
	/**
	 * this function runs RNAcofold on the Helper.tmp file and gives back the "best" hybrid as a RNAcofoldData object
	 * @return
	 */
	private static RNAcofoldData getHybrid(String file){

		List<String> duplex = RNAcofoldData.exec(GVars.RNAcofold, GVars.RNAcofoldParameter, file);
		List<RNAcofoldData> cofoldList  = null;


		cofoldList = RNAcofoldData.getCofoldDataMinBinding(duplex, 0);

		if(cofoldList.size() >= 1)
			return cofoldList.get(0);

		return null;

	}
}
