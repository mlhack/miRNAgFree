package miRNAgFreeGit;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import libs.GVars;
import libs.IO;
import libs.Read;
import libs.Sort;
import libs.Write;
import libs.Homologs;

public class GfreeHelper {

	
	
	static double fluc5p;
	static double dom2AllRatio;
	static double fluc5pMin;
	static double dom2AllRatioMin;
	static double bindings;
	static int minReadNumber; // the minimum number of reads in the 'dominant' cluster
	static int minDominantRC;
	
	static int minClusterRC=20;
	
	/**
	 * Applies the filter
	 * Obtains the microRNAs that should be used for homologous considerations
	 * @param noReads: the total number of reads in analysis
	 */
	public static void prepare(){
		
	

		
		if(Vars.mode.equals("strict")){
			IO.log(GVars.logFile, 1, "Start analysis with strict settings", true);
//			Vars.bindingPattern = readPattern(Vars.highConfPatternFile);
			setThresholdsStrict();
		}
		else{
//			Vars.bindingPattern = readPattern(Vars.allPatternFile);			
			setThresholdsLax();
		}
		bindings = Vars.bindings;
		// de-novo prediction (without homologous guide)
		if(GVars.microRNA == null){

		}
		else{
			String[] species = GVars.microRNA.split(":");

			if(new File(GVars.mature).getParent() == null)
				Read.getSpeciesMicroRNAs(GVars.libsPath+File.separator+GVars.mature,GVars.output+File.separator+"mature.fa", 
					species);
			else{
				Read.getSpeciesMicroRNAs(GVars.mature,GVars.output+File.separator+"mature.fa", 
						species);
			}
			if(new File(GVars.hairpin).getParent() == null)
				Read.getSpeciesMicroRNAs(GVars.libsPath+File.separator+GVars.hairpin,GVars.output+File.separator+"hairpin.fa", 
						species);
			else{
				Read.getSpeciesMicroRNAs(GVars.hairpin,GVars.output+File.separator+"hairpin.fa", 
						species);
			}
			if(new File(GVars.output+File.separator+"hairpin.fa").exists()){
				Vars.homologMapMature= Read.getFastaMapInverse(GVars.output+File.separator+"mature.fa");
				Vars.homologMapHairpin= Read.getFastaMapInverse(GVars.output+File.separator+"hairpin.fa");
				Vars.homologs = new Homologs(GVars.output+File.separator+"mature.fa", Vars.seedStart, Vars.seedEnd);
			}
			else{
				IO.log(GVars.logFile, 4, "Problem with microRNA file - probably no microRNAs have been found. Will quit now.", true);
				IO.warning("Problems with hairpin file - will quit now");
				System.exit(1);
			}
		}
	}
	
	/**
	 * remove last base from reads --> do not count as MM
	 * @param list
	 * @return
	 */
	public static List<String> getLengthMinusOne(List<Sort> list, int removeBases){
		
		List<String> back = new ArrayList<String>();
		for(Sort sort : list){
			back.add(sort.line.substring(0, sort.line.length()-removeBases));
		}
		return back;
	}
	
	public static void setThresholdsStrict(){
		
		if(GVars.kingdom.equals("plant")){
			fluc5pMin = 0.508474576;
			fluc5p = 0.775508096;
			

			dom2AllRatio = 0.78960353;
			dom2AllRatioMin = 0.61474819;
			
			minReadNumber = 4;
			minDominantRC = 10;
			bindings=14;

		}
		else{
			fluc5pMin = 0.52;
			fluc5p = 0.88;

			dom2AllRatioMin = 0.25;
			dom2AllRatio = 0.39;

			minReadNumber = 4;
			minDominantRC = 10;
			bindings=14;

		}
		
	}
	
	/**
	 * set the thresholds for candidate prediction
	 */
	public static void setThresholdsLax(){


		if(GVars.kingdom.equals("plant")){
			
			fluc5pMin = 0.508474576;
			fluc5p = 0.775508096;
			

			dom2AllRatio = 0.78960353;
			dom2AllRatioMin = 0.61474819;
			
			minReadNumber = 3;
			minDominantRC = 5;
			bindings=14;

		}
		else{
			
			fluc5pMin = 0.52;
			fluc5p = 0.88;

//			dom2AllRatio = 0.246292463;	
//			fluc5p = 0.6;
			dom2AllRatioMin = 0.25;
			dom2AllRatio = 0.39;
			
			minReadNumber = 3;
			minDominantRC = 5;
			bindings=14;

		}


		
	}
	
	/**
	 * read the pattern file 
	 * @param file
	 * @return
	 */
	public static Map<String,String> readPattern(String file){
		
		return  Read.readFileMap(file, 1, true, "\\t");
		
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
	//////////////////////////////
	/////////////////////////////  WRITE OUT RESULTS
	
	/**
	 * 
	 * @param resList
	 */
	public static void writeOutResults(List<GFresult> resList, int totalReads){
		
		
//		Set<Integer> detectedStar = new HashSet<Integer>();
		Set<String> detected = new HashSet<String>();
		
		int totalInLib = getTotalRCofPredicted(resList);
		String miRout = GVars.output+File.separator+"microRNAs.txt";
		String mature = GVars.output+File.separator+"mature.fa";
		String addInfo = GVars.output+File.separator+"info.txt";
		String homologLax = GVars.output+File.separator+"homologousLax.txt";
		String matureLax = GVars.output+File.separator+"maturehomologousLax.fa";
		
		boolean homo = false;

		
		new File(mature).delete();
		new File(miRout).delete();
		

		Write.writeString(miRout, 	"name\tunique reads\tread count\tcanonical read count\tRPM (lib)\tRPM (total)", false);
		Write.writeString(homologLax, 	"name\tunique reads\tread count\tcanonical read count\tRPM (lib)\tRPM (total)", false);		
		Write.writeString(addInfo, "name\thomolog miRNA\t5p_pos\t3p_pos\thomologSeq\tpattern", false);

		int c=1;
		int i = 1;
		for(GFresult res : resList){
			

	
			////
			// the id has the following order: [0]-> index of more expressed cluster (originally assumed to be 5p)
			// [1] -> the index in the cluster (-1 = dominant read, other are the index in the seq list)
			// [2] --> index of 3p cluster in clusterList 
			// [3] --> index of read in 3p Clusterobject (-1 = dominant read, other are the index in the seq list)
			String[] id = res.rnaCofoldData.name.trim().split("#");

			int index5p = Integer.parseInt(id[1]);
			int index3p = Integer.parseInt(id[3]);

			// the index in the list
//			int i5p = Integer.parseInt(id[0]);
//			int i3p = Integer.parseInt(id[2]);

//			if(detectedStar.contains(i5p) || detectedStar.contains(i3p)){
//				detectedStar.add(i5p);detectedStar.add(i3p);
//				continue;
//			}

			String name = "c_"+c;
			if(res.homologName != null){
				name = getHomoName(res.homologName);	
			}
			else{
				c++;
			}


//			detectedStar.add(i5p);detectedStar.add(i3p);
			
			if(detected.contains(name)){
				name = name+"_"+i;
				i++;
				detected.add(name);
			}
			else{
				detected.add(name);
			}
			
			String name5p = name+"-5p";
			String name3p = name+"-3p";
			
			int expr5p = res.rnaCofoldData.clust5p.getTotalRC();
			int expr3p = res.rnaCofoldData.clust3p.getTotalRC();
			if(expr5p > 2* expr3p){
				name3p = name3p+"*";
			}
			else if(expr3p > 2* expr5p){
				name5p = name5p+"*";
			}
			else{
				
			}
			
			if(res.rnaCofoldData.clust5p != null && res.rnaCofoldData.clust3p != null){
				writeSequence(mature, res.rnaCofoldData.clust5p, index5p, name5p);
				writeSequence(mature, res.rnaCofoldData.clust3p, index3p, name3p);
				
				writeExprDataExprData(miRout, res.rnaCofoldData.clust5p, name5p, totalInLib, totalReads);
				writeExprDataExprData(miRout, res.rnaCofoldData.clust3p, name3p, totalInLib, totalReads);

				res.rnaCofoldData.representDuplexStacked(res.rnaCofoldData.clust5p, res.rnaCofoldData.clust3p, 
						GVars.output+File.separator+"align"+File.separator+name+".align", Vars.overhang5p, Vars.maxReadsVisu);
//				res.rnaCofoldData.representDuplex(res.rnaCofoldData.clust5p, res.rnaCofoldData.clust3p, 
//						GVars.output+File.separator+"align"+File.separator+name+".align", Vars.overhang5p);
			}
			else{
				homo = true;
				writeSequence(matureLax, res.rnaCofoldData.clust5p, index5p, name5p);
				writeExprDataExprData(homologLax, res.rnaCofoldData.clust5p, name5p, totalInLib, totalReads);
				res.rnaCofoldData.clust5p.representCluster( GVars.output+File.separator+"align"+File.separator+name+".salign",Vars.overhang5p);
			}
			

			
			Write.writeString(addInfo, name+"\t"+res.homologName+"\t"+res.pos5p+"\t"+res.pos3p+"\t"+res.homologSequence+
					"\t"+res.knownPatternMicroRNAs, true);




		}
		
		if(!(homo)){
			new File(homologLax).delete();
//			new File(addInfo).delete();
		}
		
		
//		IO.writeToCommandLineL1("found "+c+" microRNAs");
		
	}
	
	private static void writeSequence(String file, GFcluster data, int index, String name){
		
		String seq = null;
		if(index < 0)
			seq = data.prominentSeq;
		else{
			seq = data.seqs.get(index);
		}
		Write.writeFasta(file, seq, name, "", true);

	}
	/**
	 * 
	 * @param file
	 * @param cluster
	 * @param arm
	 * @param name
	 * @param totalInLib
	 * @param totalReads
	 */
	private static void writeExprDataExprData(String file, GFcluster cluster,String name, int totalInLib, int totalReads){
	
		double rpmTotal = 1000000d*cluster.getTotalRC()/(double)totalReads;
		double rpmLib = 1000000d*cluster.getTotalRC()/(double)totalInLib;
		
		int ur = cluster.seqs.size()+1;
		Write.writeString(file, name+"\t"+ur+"\t"+cluster.getTotalRC()+"\t"+cluster.prominentRC+"\t"+rpmLib+"\t"+rpmTotal, true);

	}


	public static String getHomoName(String homoName){

		StringBuilder sb = new StringBuilder(0);
		String[] ns = homoName.split("=")[0].split("-");
		for(int i = 1; i <= 2; i++){
			if(ns.length - 1 >= i )
				sb.append(ns[i]+"-");
		}
		homoName = sb.substring(0, sb.length()-1).toString();
		return homoName;
	}
	
	/**
	 * Gives back the total number of reads in predicted microRNAs
	 * @param resList
	 * @return
	 */
	public static int getTotalRCofPredicted(List<GFresult> resList){
		
		int total = 0;
		for(GFresult res : resList){
			total += res.rnaCofoldData.clust5p.getTotalRC();
			total += res.rnaCofoldData.clust3p.getTotalRC();
		}
		return total;
	}
	
	
	
}
