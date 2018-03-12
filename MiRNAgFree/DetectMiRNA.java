package miRNAgFreeGit;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import libs.GVars;

public class DetectMiRNA {

	
	
	public static void filterCluster(List<GFcluster> clusterList){
		
		for(int i = 0; i < clusterList.size(); i++){
			double fluc = clusterList.get(i).get5pFluct() ;
			if( fluc < GfreeHelper.fluc5pMin && fluc >= 0d){
				clusterList.remove(i);
				i--;
			}
		}
	}
	
	/**
	 * This function: i) test the putative guide clusters<br> 
	 *               ii) calculate the duplex structures and retain the best as the most likely passenger sequence
	 * @param clusterList  - all read clusters as GFcluster objects
	 * @return a list of novel microRNAs (GFresult objects)
	 */
	public static List<GFresult> detectMicroRNAs(List<GFcluster> clusterList){

		
		
		List<GFresult> resList = new ArrayList<GFresult>();
		Map<Thread,TestGFcluster> threadMap = new Hashtable<Thread,TestGFcluster>();
		
		/*
		 * Go over all read clusters
		 */
		for(int i = 0; i < clusterList.size(); i++){
			
			
			if(clusterList.get(i).prominentRC < GfreeHelper.minDominantRC || 
					clusterList.get(i).getReadNumber() < GfreeHelper.minReadNumber ||
					clusterList.get(i).getTotalRC() < GfreeHelper.minClusterRC || clusterList.get(i).prominentSeq.length() < Vars.minDominantLength 
					|| clusterList.get(i).prominentSeq.length() > Vars.maxDominantLength){
				continue;
			}

			/**
			 * If the read cluster passes the criterion for 5' fluctuation and dominant to all ratio
			 */
			if(Vars.mode.equals("lax")){
				if(
						!(
								(clusterList.get(i).get5pFluct()  >= GfreeHelper.fluc5p && 
								clusterList.get(i).getDominant2All() >= GfreeHelper.dom2AllRatioMin) 

								// replace with || to make less stringent
								//					&& 
								||
								(clusterList.get(i).get5pFluct()  >= GfreeHelper.fluc5pMin && 
								clusterList.get(i).getDominant2All() >= GfreeHelper.dom2AllRatio)

								)
						){

					continue;
				}
			}
			
			
			
			if(Vars.mode.equals("strict")){
				if(
						!(
								(clusterList.get(i).get5pFluct()  >= GfreeHelper.fluc5p && 
								clusterList.get(i).getDominant2All() >= GfreeHelper.dom2AllRatioMin) 

								// replace with || to make less stringent
								&& 
								(clusterList.get(i).get5pFluct()  >= GfreeHelper.fluc5pMin && 
								clusterList.get(i).getDominant2All() >= GfreeHelper.dom2AllRatio)

								)
						){

					continue;
				}
			}
			
			/**
			 * If algorithm goes here --> guide cluster passed 
			 * --> launch a thread to calculate all possible duplex structures against all other read cluster
			 */
			
			if(threadMap.size() < GVars.p){
				TestGFcluster clust = new TestGFcluster(clusterList,i,GVars.tmp+"_"+i);
				Thread t = new Thread(clust);
				t.start();
				threadMap.put(t, clust);
			}
			else{
				checkThreads(threadMap, resList);
				i--;
			}
		}
		
		// wait until last thread has finished
		while(threadMap.size() > 0){
			checkThreads(threadMap, resList);
		}

		return resList;
	}

	public static void checkThreads(Map<Thread,TestGFcluster> threadMap, List<GFresult> resultList){

		for(Iterator<Map.Entry<Thread,TestGFcluster>> it = threadMap.entrySet().iterator(); it.hasNext();){
			Map.Entry<Thread,TestGFcluster> entry = it.next();
			if(entry.getKey().isAlive()){

			}
			else{
				if(entry.getValue().result != null){
					resultList.add(entry.getValue().result);
				}

				it.remove();
			}
		}
	}

}
