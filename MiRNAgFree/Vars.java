package miRNAgFreeGit;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import libs.Classification;
import libs.Exec;
import libs.GVars;
import libs.IO;
import libs.Read;
import libs.Util;
import libs.Homologs;

public class Vars {

	
	static public int maxParamterLen = 30;
	static public int maxWindowLen = 70;
	static public String version = "04/17";
	static public Set<String> ap = new HashSet<String>(); // all parameters are stored in this Hash --> check the input
	static public String alignOut;  // the output for the alignment files
	static Map<String,String> bindingPattern; // map for the hybrid pattern 
	static public String allPatternFile; // pattern of all microRNAs
	static public String highConfPatternFile; // pattern of only high-confidence pattern 
	
	static public String modelFile = null; // the model file used for classification
	static public String modelFileAnimal = "rf_animal.model"; // the model file used for classification
	static public String modelFilePlant = "rf_plant.model"; // the model file used for classification

	static public Classification classification = null;
	
	static public List<String> libsFilter; // the names of the libraries that should be used for filtering
	static public Map<String,String> homologMapMature; // the map that holds the putative homologous mature sequences
	static public Map<String,String> homologMapHairpin; // themap that holds the putative homologous hairpin sequence
	static public Homologs homologs =null;
	
	static public int seedStart = 2; // start position of the seed
	static public int seedEnd = 8; // end position of the seed
	static public String mode = "strict"; // the mode: strict; lax (and others that might be implemented in the future
	static public String procMode = "lax"; // the processing; strict = a strict 2nt overhang is requiered; lax = a 1 nt tolerance is allowed
	
	static public int minDominantLength = 20; // default minimum length of dominant read
	static public int maxDominantLength = 23; // default maximum length of dominant read
	
	static public int overhang3p = 3; // the 3p overhang for the read clustering
	static public int overhang5p = 3; // the 5p overhang for the read clustering
	static public int guideMM = 2; // max. number of mismatches against the guide strand
	static public int passMM = 6; // max. number of mismatches against the passenger strand
	
	static int nrNonDomGuide = -1	; //the maximum number of non-dominant reads that should be checked for hybrids in the dominant cluster (guide)
	static int nrNonDomStar = -1; // number of putative passenger strands that are checked (the passenger sequence must be within the first 'nrNonDomStar' reads in the passenger cluster

	static int maxReadsVisu = 10; // the maximum number of reads for visual output
	
	static int patternLength = 5; // 
	static double scoreThreshold = 0.9; // the probability that the instance is positive
	
	static int mm = 1;
	static int bindings = 14;
	
	static boolean perfectHybrid = true; // if perfect hybrids should be accepted

	/**
	 * 1) reads the parameters either from the command line or from a config file
	 * 2) sets the default values for those that do not appear
	 * 3) prepare the output directories
	 * @param args: the command line 
	 */
	public static void getParameters(String args[]){
		
		/**
		 * Define existing parameters
		 */
		GVars.setExistingParameters();
		setExistingParameters();
		ap.addAll(GVars.ap);
		
		
		
		Map<String,List<String>> info = null;
		////// SET THE PARAMETERS  //////////////
		if(args.length > 0 && new File(args[0]).exists()){
			// Read configuration script
			info = IO.readConfigFile(args[0]);
		}
		else{
			info = getCommandLineParameter(args);
		}
		
		// set general default values
		GVars.addDefaults(info);
		// set sRNAbench specific default values
		addDefaults(info);
		
		// prepare the output 
		new File(GVars.stat).mkdir();
		new File(Vars.alignOut).mkdir();		
		if(GVars.graphics)
			new File(GVars.graphs).mkdir();
		
		sanityCheck();
		
	}
	
	public static void sanityCheck(){
		
		String inst = Exec.checkIfProgram(GVars.RNAcofold);
		if(inst == null)
		{
			System.out.println(GVars.RNAcofold+" was not found to be installed. Please install first. Will quit now!");
			System.exit(1);
		}
		else{
			System.out.println("Found "+GVars.RNAcofold+" in "+inst);
		}
//		System.out.println("implement sanity check in Vars");
	}
	/**
	 * 
	 * @param info
	 */
	public static void addDefaults(Map<String,List<String>> info){
		
		
		GVars.stat = GVars.output+File.separator+"stat";
		alignOut = GVars.output+File.separator+"align";
		
//		allPatternFile = GVars.dbPath + File.separator + "pattern" + File.separator + "pattern_all.txt";
//		highConfPatternFile =  GVars.dbPath + File.separator + "pattern" + File.separator + "pattern_highConf.txt";

		
		if(info.containsKey("libsFilter")){
			libsFilter = new ArrayList<String>();
			for(String l : info.get("libsFilter"))
				libsFilter.add(l);
		}	
		if(info.containsKey("mode")){
			mode = info.get("mode").get(0);
		}
		
		if(info.containsKey("overhang3p")){
			overhang3p = Integer.parseInt(info.get("overhang3p").get(0));
		}
		if(info.containsKey("overhang5p")){
			overhang5p = Integer.parseInt(info.get("overhang5p").get(0));
		}
		
		
		if(info.containsKey("minReadLength")){
			GVars.minReadLength = Integer.parseInt(info.get("minReadLength").get(0));
		}
		else{
			IO.writeToCommandLineL1("Set min. read length to 19!");
			GVars.minReadLength = 19;
		}
		if(info.containsKey("maxReadLength")){
			GVars.maxReadLength = Integer.parseInt(info.get("maxReadLength").get(0));
		}
		else{
			IO.writeToCommandLineL1("Set max. read length to 25!");
			GVars.maxReadLength = 25;
		}
		
		if(info.containsKey("minDominantLength")){
			Vars.minDominantLength = Integer.parseInt(info.get("minDominantLength").get(0));
		}
		if(info.containsKey("maxDominantLength")){
			Vars.maxDominantLength = Integer.parseInt(info.get("maxDominantLength").get(0));
		}

		
		if(info.containsKey("guideMM")){
			guideMM = Integer.parseInt(info.get("guideMM").get(0));
		}
		if(info.containsKey("passMM")){
			passMM = Integer.parseInt(info.get("passMM").get(0));
		}

		if(info.containsKey("nrNonDomGuide")){
			nrNonDomGuide = Integer.parseInt(info.get("nrNonDomGuide").get(0));
			nrNonDomGuide = nrNonDomGuide - 1;
		}
		if(info.containsKey("nrNonDomStar")){
			nrNonDomStar = Integer.parseInt(info.get("nrNonDomStar").get(0));
			nrNonDomStar = nrNonDomStar - 1;
		}

		if(info.containsKey("mm")){
			mm = Integer.parseInt(info.get("mm").get(0));
		}

		if(info.containsKey("procMode")){
			procMode = info.get("procMode").get(0);
		}
		
		if(info.containsKey("maxReadsVisu")){
			maxReadsVisu = Integer.parseInt(info.get("maxReadsVisu").get(0));
		}
		if(info.containsKey("scoreThreshold")){
			scoreThreshold = Double.parseDouble(info.get("scoreThreshold").get(0));
		}
		
		
		if(info.containsKey("bindings")){
			bindings = Integer.parseInt(info.get("bindings").get(0));
		}

		
		if(info.containsKey("modelFile")){
			modelFile = info.get("modelFile").get(0);
		}
		else if(GVars.kingdom.equals("animal")){
			IO.writeToCommandLineL1("Will use animal model for prediction. For plants: kingdom=plant");
			modelFile = GVars.dbPath + File.separator + "pattern" + File.separator + modelFileAnimal;
		}
		else if(GVars.kingdom.equals("plant")){
			modelFile = GVars.dbPath + File.separator + "pattern" + File.separator + modelFilePlant;	
			
			perfectHybrid = false;
			if(info.containsKey("maxReadLength")){

			}
			else{
				GVars.maxReadLength = 23;
			}

		}
		else{
			IO.warning("Fatal error for parameter kingdom which must be either kingdom=animal or kingdom=plant !");
			System.exit(1);
		}
//		classification = new Classification(modelFile);

		if(info.containsKey(perfectHybrid)){
			perfectHybrid = Boolean.parseBoolean(info.get("perfectHybrid").get(0));			
		}
		
		
		if(GVars.microRNA != null){
			
			String hairpinFile = GVars.libsPath+File.separator+GVars.hairpin;
			if(new File(GVars.hairpin).exists()){
				hairpinFile = GVars.hairpin;
			}
			
			int parsed = Read.getSpeciesMicroRNAs(hairpinFile,GVars.output+File.separator+"hairpin.fa", 
					GVars.microRNA.split(":"));
			Vars.homologMapHairpin = Read.getFastaMap(GVars.output+File.separator+"hairpin.fa");
			
			String matureFile = GVars.libsPath+File.separator+GVars.mature;
			if(new File(GVars.mature).exists()){
				matureFile = GVars.mature;
			}
			
			int parsedM = Read.getSpeciesMicroRNAs(matureFile,GVars.output+File.separator+"mature.fa", 
					GVars.microRNA.split(":"));
			Vars.homologMapMature = Read.getFastaMap(GVars.output+File.separator+"mature.fa");

			
			Vars.homologs = new Homologs(GVars.output+File.separator+"mature.fa",seedStart,seedEnd);
		}

	}
	/**
	 * Add the existing parameters to the 'ad' Set
	 */
	public static void setExistingParameters(){
		
		ap.add("alignOut");
		ap.add("mode");
		ap.add("libsFilter");
		ap.add("overhang3p");
		ap.add("overhang5p");
		
		

		ap.add("bindings");

		ap.add("guideMM");
		ap.add("passMM");
		ap.add("nrNonDomGuide");
		ap.add("nrNonDomStar");

		ap.add("procMode");

		ap.add("maxReadsVisu");

		ap.add("modelFile");
		
		ap.add("scoreThreshold");
		ap.add("mm");
		
		ap.add("perfectHybrid");

		ap.add("minDominantLength");
		ap.add("maxDominantLength");
		
	}
		
	/**
	 * reads command line and writes it into the info Map in Helper
	 * @param args
	 * @return
	 */
	public static Map<String,List<String>> getCommandLineParameter(String[] args){
		
		Map<String,List<String>> back = new Hashtable<String,List<String>>();
		if(args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			showOptions();
		}
		else if(args[0].equals("-v")){
			IO.warning("You are using version: "+version);
			System.exit(1);
		}
		for(String arg : args){
//			System.out.println(arg);
			String[] f = arg.split("=");
			if(f.length == 2){
				if(ap.contains(f[0])){
					Util.setInfo(back,f[0], f[1]);
				}
				else{
					IO.warning("The following parameter given on the command line does not exist: "+f[0]);
					System.exit(1);
				}
			}
			else{
				IO.warning("Fatal error in parameters: "+arg+"");
				IO.warning("The parameters must be in parameterName=value format");
				System.exit(1);
			}
		}
		
		return back;
	}

	/**
	 * Show the most important input options
	 */
	public static void showOptions(){

		IO.writeToCommandLineBlockOpen("GENERAL THINGS");
		IO.writeToCommandLineL1("The parameters must be given in this format: parameter=value");
		IO.writeToCommandLineL2( "For example input=sample.fastq");
		
		IO.writeToCommandLineBlockOpen("");

		System.out.println("Some basic parameters are explained below. Please see the manual for more details");
		System.out.println("------------------------\n");

		Util.printArguments("input=<file>","The path to the input file (supported formats: fastq, read/count, fasta). This is the only mandatory parameter.", maxParamterLen,maxWindowLen);
		Util.printArguments("output=<folder>","The name of the oputput folder. Default: output=/opt/sRNAtoolboxDB/out ", maxParamterLen,maxWindowLen);
		Util.printArguments("dbPath=<path>","Path to sRNAtoolbox database. Default: dbPath=/opt/sRNAtoolboxDB ",maxParamterLen,maxWindowLen);		

		Util.printArguments("libsFilter=<path>","Path to libraries that should be used for reads filtering. ",maxParamterLen,maxWindowLen);		

		
		
		Util.printArguments("microRNA=<String>","Specify the microRNAs that should be used as a set of homologous microRNA. The short name of species (like used in miRBase). For example" +
				" microRNA=hsa or microRNA=mmu are used. Several short species names can be given like microRNA=hsa:mmu:rno:ebv", maxParamterLen,maxWindowLen);

		Util.printArguments("mode=<strict,lax>","mode=strict: Both, for the thresholds and the patterns the high confidence parameters are used. Default mode=lax ",maxParamterLen,maxWindowLen);
		Util.printArguments("procMode=<strict,lax>","mode=strict: a strict 2nt overhang is requiered; lax = a 1 nt tolerance is allowed Default: procMode=lax",maxParamterLen,maxWindowLen);

		Util.printArguments("mm=<int>","the number of mismatches within the microRNA cluster (stack). Default: mm=1",maxParamterLen,maxWindowLen);

		Util.printArguments("bindings=<int>","the minimum number of bindings that a hybrid must have Default: bindings=14. ",maxParamterLen,maxWindowLen);

		Util.printArguments("nrNonDomGuide=<int>","the number of non-dominant "
				+ "reads that should be checked for hybrids in the dominant cluster (guide). "
				+ "Default: nrNonDomGuide=-1 (only the dominant read of the most expressed cluster is considered) ",maxParamterLen,maxWindowLen);
		
		Util.printArguments("nrNonDomStar=<int>","the number of non-dominant "
				+ "reads that should be checked for hybrids of the star cluster (passanger). "
				+ "Default: nrNonDomStar=-1 (only the dominant read of the passenger cluster is considered) ",maxParamterLen,maxWindowLen);

//		Util.printArguments("scoreThreshold=<double>","the probability that the analysed hybrid corresponds to a microRNA hybrid. Default: scoreThreshold=0.9",maxParamterLen,maxWindowLen);

		
		Util.printArguments("kingdom=<String>","Indicate if the data is from animal (kingdom=animal) or plant (kingdom=plant). This paramter "
				+ "affects only the prediction of novel microRNAs. Default: kingdom=animal", maxParamterLen,maxWindowLen);

		Util.printArguments("perfectHybrid=<true,false>","true --> the method allows perfect (no internal loops or bulges) hybrids. By default, perfectHybrid=true. If kingdom=plant, "
				+ "this parameter is set to false by default to avoid prediction NAT-siRNAs, however, it can be overwritten by the user on the command line. ", maxParamterLen,maxWindowLen);
		
		Util.printArguments("minDominantLength=<int>","The minimum length the dominant sequence (canonical sequence) must have. "
				+ ""
				+ "Default: minDominantLength=20 ",maxParamterLen,maxWindowLen);

		Util.printArguments("maxDominantLength=<int>","The maximum length the dominant sequence (canonical sequence) can have. "
				+ ""
				+ "Default: minDominantLength=23 ",maxParamterLen,maxWindowLen);

		



		System.exit(1);

	}

	/**
	 * Print the sRNAbench  header
	 */
	public static void welcome(){


		System.out.println("\n\n"+Util.getCharString(IO.outWidth, '*'));
		System.out.println(Util.getCharString(IO.outWidth, '*'));
		System.out.println(Util.getCharString(10, '*')+Util.getVoidString(IO.outWidth-20)+Util.getCharString(10, '*')+"\n");

		System.out.println("           miRNAgFree version "+version);
		System.out.println("           Computational Epigenomics Group ");
		System.out.println("           Genetics Department, University of Granada, Spain  ");
		System.out.println("           For more information, please visit:  http://bioinfo2.ugr.es \n");


		System.out.println(Util.getCharString(10, '*')+Util.getVoidString(IO.outWidth-20)+Util.getCharString(10, '*'));
		System.out.println(Util.getCharString(IO.outWidth, '*'));
		System.out.println(Util.getCharString(IO.outWidth, '*')+"\n");

	}
	


}
