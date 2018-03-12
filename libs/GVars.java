package libs;

import java.io.File;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;



/**
 * General Variables for sRNAtoolbox tools
 * @author Michael
 *
 */
public class GVars {

	
	

	
	// existing parameters
	static public Set<String> ap = new HashSet<String>(); // all parameters are stored in this Hash --> check the input
	static public Map<String,List<MapData>> genomeMap= new Hashtable<String,List<MapData>>(); // store the mappings to the genome
//	static public Map<Integer,List<MapData>> readIDmap= new Hashtable<Integer,List<MapData>>(); // store the mappings to the genome
	/////////////////////////////////////////////////////////////////////////////////
	// other programs 
	static public String aligner = "bowtie"; 
	static public String bowtie = "bowtie";
	static public String bowtieBuild = "bowtie-build";
	static public String bowtieInspect = "bowtie-inspect";
	static public String RNAfold = "RNAfold";
	static public String RNAduplex = "RNAduplex";
	static public String RNAcofold = "RNAcofold";
	static public String RNAfoldParameter = " --noPS  ";
	static public String RNAcofoldParameter = "  --noPS  ";
	static public String RNAduplexParameter = "  --noLP ";
	static public String TargetFinder = "perl /usr/local/bin/targetfinder_mod.pl";
	

	///////////////////////////////////////////////////////////////////////////////////
	//// general IO paths
	static public String output; // the output folder
	static public String input;  // the input file 
	static public String dbPath; // the path to the database
	
	/// parameters needed in preprocessing and/or parameters that influence the behaviour
	static public String species; // the name of the bowtie index
	
	static public int p = 4; // the number of cpus asigned
	
	static public String sep = "#"; // the separator used in fasta file input

	/////////////////////////////////////////////////////////////////////
	/// LIBRARIES
	static public String microRNA = null; // the microRNA string microRNA=hsa:xxx
	static public String homolog = null; // the microRNA string microRNA=hsa:xxx
	static public String hairpin = "hairpin.fa";
	static public String mature = "mature.fa";	
	
	/////////////////////////////////////////////////////////////////////
	/// SOLID 
	static public boolean solid = false; // the input is SOLID
	// these two are only internal 
	static public String colorIndex = ""; // void for nucleotide space, "_c" for color space
	static public String colorFlag = ""; // the parameter that needs to be given, i.e. -C (for bowtie)

	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	///// internal names  - those cannot be given by command line parameters
	
	static public String genomeZip; // the file name of the genome mapping file
	static public boolean doAlignment = true; // i.e. the input is not from previously aligned data
	static public String inputType = null; // will be overwritten by Preprocessing functions
	// INTERNAL PARAMETERS OR PATHS  -- database related
	static public String libsPath; // the path to the libs folder
	static public String seqOBJ; // the path to the seqOBJ folder
	static public String index; // the path to the bowtie index
	static public String origInput; // the original input file - this cannot be set by means of parameters
	static public String tmp; // the temporal file 
	static public String logFile; // the new unified log file name
	static public String parameters; // the parameters file
	////// output folders 
	static public String stat; // the path to the 'stat' folder
	static public String graphs; // the path to the graphs folder 
	static public String bigWig; // the output folder for bedGraph and bigWig files
	//////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////
	// quality parameters
	
	static public int phred = 33; // the phred score encoding --> this number will be rested 
	static public int minQ = 20; // the worst base call quality --> if worse that that, the read will be trimmed at this base
	static public String qualityType = null; // the type of quality control: min=all bases must be above the minimum;
	static public double maxQfailure = 0;
	 // mean=Q must be above this mean 
	// adapter parameters
	static public int adapterStart = 0; // where to start the search for the adapter
	static public int adapterMM = 1; // number of mismatches allowed in the first 'adapterMinLength' bases
	static public int adapterMinLength = 10; // minimum alignment length between adapter sequence and read
	static public boolean holdNonAdapter = false; // if set to true --> untrimmed reads will be hold for the analysis
	static public boolean recursiveAdapterTrimming = false; // if true --> do an recursive search for the adapter sequence (making adapterMinLength smaller in each round)
	static public boolean writeNonAdapter = false; // write out a file with the untrimmed reads (default FALSE)
	static public String adapter; // the adapter sequence
	static public boolean adapterTrimmed = true;
	static public boolean guessAdapter = false;
	
	//////////////////////////////////////////////////////////////////////////////////
	// Barcodes and random adapters
	static public int removeBarcode = 0; // the number of bases at the 5' end of the read that will be removed
	static public int remove3pBases = 0; // remove X nucleotides after adapter trimming 

	///////////////////////////////////////////////////////////////////////////////////////////
	/// length and count thresholds
	static public int minReadLength = 15; // the minimum read length --> shorter reads are filtered out
	static public int maxReadLength = 5000; // the maximum read length
	static public int minRC = 1; // the minimum read count
	
	/////////////////////////////////////////////////////////////////////////////////
	// MAPPING PARAMETERS
	static public int mm = 1; // number of mismatches  - applies to both alignment types
	static public int seed = 19; // the seed if seed alignment is activated
	static public String alignType = "n"; // the alignment type - can be either 'n' (Bowtie seed alignment) or 'v' (full read length alignment)
	static public String bowtieReportType = "-a -m";
	static public String bowtieReportCount = 10+"";
	static public String libAlignType = "v"; // the alignment type that should be used for the annotations / libraries (for example to map hairpins to the genome 
	static public String bowtieAdd = "--best --strata";
	static public int chunkmbs = 256; // bowtie parameter

	static public boolean direct = true;
	
	static public int libSeed = 20; // the seed used to map the libraries in fasta format 
	static public String mappingOrientation = ""; // this is an internal variable which will be passed to bowtie

	// the following 3 variables overwrite mappingOrientation in the different analysis steps
	static public String microRNAmappingOrientation = "";// --norc or --nofw --> if empty then both /
	static public String libsmappingOrientation = "";// --norc or --nofw --> if empty then both /
	static public String tRNAmappingOrientation = "";// --norc or --nofw --> if empty then both /


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//// expression profiling parameters
	static public int winUpMir = 3; // upstream window for 1) isomiRs andd 2) prediction of novel microRNAs
	static public int winDownMir = 5; // downstream window

	static public int winUpLibs = 0; // upstream window of libraries 
	static public int winDownLibs = 0; // downstream window
	
	static public int winUpPre = 0; // upstream window of libraries 
	static public int winDownPre = 0; // downstream window
	
	static public boolean hierarchical = true; // hierarchical detection --> remove after profiling one library

	//////////////////////////////////////////////////////////////////////////////////////////
	//// 	ISOMIR PARAMETERS
	static public int isomiRseed = 18; // the seed to determine the isomiRs
	static public boolean isoMiR = false; // calculate isomiRs
	static public boolean isoLibs = false; // calculate isoElement for others than microRNAs
	static public int minRCiso = 10;  // minimum RC for isomiR profiling (if RC is smaller, this RNA is not considered)
	static public boolean fullIsoStat = false;

	//////////////////////////////////////////////////////////////////
	///// OUTPUT OPTIONS
	static public boolean graphics = false; // make graphics
	static public boolean plotLibs = false; // write out alignments and sec. structures for others that microRNAs
	static public boolean plotMiR = true; // write out miR alignments
	static public boolean plotSec = false; // write out miR alignments

	static public boolean bedGraph = false; // make bedGraph files of genome mapping
	static public String bedGraphMode = "FA"; // FA (full assignment) or MA (multiple assignment)
	static public String bedGraphIntervals = null; // the length intervals that should be used

	static public int minRCplotLibs = 200; // min. RC to write out alignment file
	static public int minRCplotMiR = 20;  // min. RC to write out alignment file
	
	static public int maxLenSecStruc = 200; // maximum length of RNA to calculate the secondary structure
	static public int maxLenPlot = 2000; // the maximum length that a reference sequence can have

	
	static public String sRNAassoc = null; // use a association file
	static public String libsStringTypes = null; //to make a specific summary file
	static public String libsStringNames = null; 
	
	static public boolean remove = true; // remove some of the result files
	static public boolean removeMicroReads = true; // remove the microRNA reads after mapping


	///////////////////////////
	/// analysis types 
	
	static public String fileType = "-f"; // the input file type for bowtie mapping 
	
	
	static public boolean RNAseq = false; // simple RNAseq analysis
	static public boolean noMatureAnnot = false; // do not use mature sequences for annotations (most expressed read will be used as mature guide)
	static public boolean calcFeaturesMiR = false; // calculates the processing and secondary structure features of the microRNAs
	
	static public boolean calcFeaturesDominant = false; // takes the most expressed read for feature



	static public int iterative5pTrimming = -1; // interative 5' read trimming 
	//////////////////////////////////////////////////////////////////



	/// the temporal mapping parameters
	static public int tempMM = 1;
	static public int tempSeed = 19;
	static public String tempAlingType = "n";
	static public String tempBowtieReportType = "-a -m"; 
	static public String tempBowtieReportCount = 10+"";
	

	static public String rnaSeqLib; // the transcripts for RNA-seq profiling

	

	///// MAKE GENOME DISTRIBUTION
	static public boolean writeGenomeDist = true; // make genome distribution
	static public boolean splitToSpecies = true; // split the genome mappings to the different species (write out the files)
	static public String assocFile; // a file that groups the different species into a higher order classification
	static public boolean chromosomeLevel = false; // make distributions at a chromosome level
	static public String mainSpecies = null; // 
	static public boolean genomeDistunique = false; // make genome distribution only for those that are uniquely mapped to one taxonomic order
	
	static public boolean indChr = false; // use the name of the individual chromosome or the species?
	static public boolean chrMappingByLength = false; // only applies if writeGenomeDist is set to TRUE; writes out the chromosome distribution as a function of read length
	static public String chromosomes = null; // the chromosomes that should be used for 'chrMappingByLength'
	
	static public boolean chrLength = false;
	

	
	/////////////////////////////////////////////////////////////////////////////////////
	/// library parameters
	static public int base = 1; // for bed file input: 1 --> add 1 to start (convert UCSC 0-based to 1-based); 0 --> do nothing
	static public String gtfAtribute = "transcript_id"; // either 
	static public String gtfFeature = "exon";
	
	static public int matureMM = 0; // the number of mismatches of the reference sequence to the genome
	static public int hairpinMM = 1;	// the number of mismatches of the reference sequence to the genome
	static public int matureHomologMM = 2; // number of mismatches when mapping homologous microRNAs to the genome
	static public int libsMM = 1;  // number of mismatches when mapping the libs to the genome
	
	//////////////////////////////////////////////////////////////////7
	///  CLUSTER DETECTION  
	static public int clusterUp = 4 ;
	static public int clusterDown = 4;
	static public int minRCcluster = 1;
	static public int minReadsInCluster=3;
	

	
	
	static public boolean nonRed = false;
	
	
	static public int foldBackTol = 0; // the foldback tolerance (for example 451a folds back in 1nt onto itself)
	
	
	//////////////////////////////////////////////////////////////7
	/////  
	////  	prediction parameters
	static public String kingdom = "animal";
	static public boolean predict = false;
	static public int flankNovel = 11; // the flanking region added to the pre-miRNA
	static public int maxDistNovel = 0; // the maximal distance between the two read stacks that might correspond to the novel microRNAs 
	static public int maxDistAnimal = 60;
	static public int maxDistPlant = 180;
	static public String novelName = "new";
	static public int startCountNovel = 1;
	static public String novelHomolog = null;
	
	static public int maxNonDominant = 2;
	static public int maxNonDominantStar =5;
	
	
	static public int seedStart = 2;
	static public int seedEnd = 8;
	static public int seedShiftSize = 0;
	static public int seedMM = 0;
	
	//////////////////////////////////////////////////////////////////////////////////////////

	public static void setExistingParameters(){
		
		
		/////////////////////////////////////////////////////////////////////////////////
		// other programs 
		ap.add("aligner");
		ap.add("bowtie");
		ap.add("bowtieBuild");
		ap.add("inspect");
		ap.add("RNAfold");
		ap.add("RNAduplex");
		ap.add("RNAcofold");
		ap.add("RNAfoldParameter");
		ap.add("RNAcofoldParameter");
		ap.add("RNAduplexParameter");
		ap.add("TargetFinder");
		
		///////////////////////////////////////////////////////////////////////////////////
		//// general IO paths
		ap.add("output");
		ap.add("input");
		ap.add("tmp");
		ap.add("logFile");
		ap.add("parameters");
		
		///////////////////////////////////////////////////////////////////////////////////////
		// INTERNAL PARAMETERS OR PATHS  -- database related
		ap.add("libsPath");
		ap.add("seqOBJ");
		ap.add("index");
		ap.add("dbPath");
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		/// parameters needed in preprocessing and/or parameters that influence the behaviour
		ap.add("species");
		ap.add("doAlignment");
		ap.add("inputType");
		ap.add("RNAseq");		
		ap.add("noMatureAnnot");		
		ap.add("calcFeaturesMiR");		
		ap.add("calcFeaturesDominant");
		ap.add("iterative5pTrimming");
		
		//////////////////////////////////////////////////////////////////////////////////
		// quality parameters
		// adapter parameters
		ap.add("phred");
		ap.add("minQ");
		ap.add("qualityType");
		ap.add("maxQfailure");

		//////////////////////////////////////////////////////////////////////////////////
		// adapter parameters
		ap.add("adapterStart");
		ap.add("adapterMM");
		ap.add("adapterMinLength");
		ap.add("holdNonAdapter");
		ap.add("recursiveAdapterTrimming");
		ap.add("writeNonAdapter");
		ap.add("adapter");
//		ap.add("adapterTrimmed");
		ap.add("guessAdapter");

		/////////////////////////////////////////////////////////////////////////////////
		// MAPPING PARAMETERS
		ap.add("noMM");
		ap.add("mm");
		ap.add("seed");
		ap.add("alignType");
		ap.add("bowtieReportType");
		ap.add("bowtieReportCount");
		ap.add("libAlignType");
		ap.add("bowtieAdd");
		ap.add("direct");
		
		ap.add("chunkmbs");
		ap.add("p");
		ap.add("libSeed");
		ap.add("solid");
		ap.add("mappingOrientation");
		ap.add("microRNAmappingOrientation");		
		ap.add("libsmappingOrientation");		
		ap.add("tRNAmappingOrientation");		
		
		ap.add("fileType");
		
		/////////////////////////////////////////////////////////////////////
		/// LIBRARIES
		ap.add("microRNA");
		ap.add("homolog");
		ap.add("hairpin");
		ap.add("mature");
		ap.add("rnaSeqLib");
		
		
		
		
		///////////////////////////////////////////////////////////////////////////////////////////
		/// preprocessing 
		ap.add("minReadLength");
		ap.add("maxReadLength");
		ap.add("minRC");
		ap.add("removeBarcode");
		ap.add("remove3pBases");
		ap.add("sep");

		
		//////////////////////////////////////////////////////////////////////////////////////////
		//// 	ISOMIR PARAMETERS
		ap.add("isomiRseed");
		ap.add("isoMiR");
		ap.add("isoLibs");
		ap.add("minRCiso");

		ap.add("fullIsoStat");
		
		//////////////////////////////////////////////////////////////7
		/////  
		////  	prediction parameters
		ap.add("prediction");
		ap.add("kingdom");
		ap.add("predict");
		
		ap.add("flankNovel");
		ap.add("maxDistNovel");
		ap.add("novelName");
		ap.add("startCountNovel");		
		ap.add("novelHomolog");
		//////////////////////////////////////////////////////////////////
		///// OUTPUT OPTIONS
		ap.add("graphics");
		ap.add("plotLibs");
		ap.add("plotMiR");
		ap.add("plotSec");		
		
		ap.add("bedGraph");
		ap.add("bedGraphMode");
		ap.add("bedGraphIntervals");
		ap.add("minRCplotLibs");
		ap.add("minRCplotMiR");
		ap.add("maxLenSecStruc");
		ap.add("maxLenPlot");

		ap.add("sRNAassoc");
		ap.add("libsStringTypes");
		ap.add("libsStringNames");

		
		
		//////////////////////////////////////////////////////////////////
		///// MAKE GENOME DISTRIBUTION
		
		ap.add("writeGenomeDist");
		ap.add("splitToSpecies");
		ap.add("assocFile");

		ap.add("mainSpecies");
		ap.add("chromosomeLevel");
		ap.add("genomeDistunique");
		ap.add("indChr");
		ap.add("chrMappingByLength");
		ap.add("chromosomes");
		ap.add("chrLength");
		
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//// expression profiling parameters
		ap.add("winUpMir");
		ap.add("winDownMir");
		ap.add("winUpLibs");
		ap.add("winDownLibs");
		
		ap.add("winUpPre");
		ap.add("winDownPre");
		ap.add("hierarchical");
		



	}
	
	

	/**
	 * add the default parameters
	 * @param info
	 */
	public static void addDefaults(Map<String,List<String>> info ){

		/////////////////////////////////////////////////////////////////////////////////
		// other programs 
		if(info.containsKey("aligner")){
			aligner = info.get("aligner").get(0);
		}
		if(info.containsKey("bowtie")){
			bowtie = info.get("bowtie").get(0);
		}
		if(info.containsKey("RNAfold")){
			RNAfold = info.get("RNAfold").get(0);
		}
		if(info.containsKey("RNAduplex")){
			RNAduplex = info.get("RNAduplex").get(0);
		}
		if(info.containsKey("RNAcofold")){
			RNAcofold = info.get("RNAcofold").get(0);
		}

		////////////////////////////////////////////////
		///////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////////////////////////////
		// INTERNAL PARAMETERS OR PATHS  -- database related
		
		if(!(info.containsKey("dbPath"))){
			dbPath = "/opt/sRNAtoolboxDB";
			libsPath = dbPath+File.separator+"libs";
			seqOBJ = dbPath+File.separator+"seqOBJ";
			index = dbPath+File.separator+"index";
			IO.writeToCommandLineL1("Using the default database path");
//			IO.warning("Using the default Please specify the database path: dbPath=path ");
//			System.exit(0);
		}
		else{
			dbPath = info.get("dbPath").get(0);
			libsPath = dbPath+File.separator+"libs";
			seqOBJ = dbPath+File.separator+"seqOBJ";
			index = dbPath+File.separator+"index";
		}
		if(!(info.containsKey("input"))){
			IO.warning("Please specify the input file: input=file ");
			System.exit(0);
		}
		else{
			input = info.get("input").get(0);
			origInput = input;
		}
		
		///////////////////////////////////////////////////////////////////////////////////
		//// general IO paths
		if(!(info.containsKey("output"))){
			IO.warning("Using default output!");
			File f = new File(input);
			String base = f.getName().split("\\.")[0];
			output = dbPath+File.separator+"out"+File.separator+base;
		}
		else{
			output = info.get("output").get(0);	
		}
		tmp = output+File.separator+"tmp";
		logFile = output+File.separator+"logFile.txt";
		parameters = output+File.separator+"parameters.txt";	
		File o = new File(GVars.output);
		if(!(o.isDirectory())){
			boolean mkdDircheck = o.mkdirs();
			if(!(mkdDircheck)){ 
				IO.warning("Output directory "+GVars.output+" cannot be created! Will quit now!");
				System.exit(1);
			}
		}
		new File(logFile).delete();

		/////////////////////////////////////////////////////////////////////////////////////////////
		/// parameters needed in preprocessing and/or parameters that influence the behaviour

		if(info.containsKey("species")){
			species = info.get("species").get(0);
		}
		else{
			IO.log(GVars.logFile, 1, "No species is given, set bowtie parameters by default to '-a'", true); 
			GVars.bowtieReportType = "-a";
			GVars.bowtieReportCount = "";

		}
		if(info.containsKey("doAlignment")){
			doAlignment = Boolean.parseBoolean(info.get("doAlignment").get(0));
		}
		if(info.containsKey("inputType")){
			inputType = info.get("inputType").get(0);
		}
		
		if(info.containsKey("RNAseq")){
			RNAseq = Boolean.parseBoolean(info.get("RNAseq").get(0));
		}
		if(info.containsKey("noMatureAnnot")){
			noMatureAnnot = Boolean.parseBoolean(info.get("noMatureAnnot").get(0));
		}		
		

		//////////////////////////////////////////////////////////////////////////////////
		// quality parameters

		if(info.containsKey("minQ")){
			minQ = Integer.parseInt(info.get("minQ").get(0));
		}
		if(info.containsKey("phred")){
			phred = Integer.parseInt(info.get("phred").get(0));
		}
		if(info.containsKey("maxQfailure")){
			maxQfailure = Integer.parseInt(info.get("maxQfailure").get(0));
		}
		if(info.containsKey("qualityType")){
			qualityType = info.get("qualityType").get(0);
			if(qualityType.equals("min") || qualityType.equals("mean")){
				
			}
			else{
				IO.warning("qualityType must be either min (eliminate reads if X bases - default X=0 - are below the minimum "
						+ "read quality) or mean (the mean quality needs to be above the minimum quality). "
						+ "\nWill quit now! ");
				System.exit(1);
			}
		}



		//////////////////////////////////////////////////////////////////////////////////
		// adapter parameters
		if(info.containsKey("adapterStart")){
			adapterStart = Integer.parseInt(info.get("adapterStart").get(0));
		}
		if(info.containsKey("adapterMM")){
			adapterMM = Integer.parseInt(info.get("adapterMM").get(0));
		}
		if(info.containsKey("adapterMinLength")){
			adapterMinLength = Integer.parseInt(info.get("adapterMinLength").get(0));
		}
		if(info.containsKey("holdNonAdapter")){
			holdNonAdapter = Boolean.parseBoolean(info.get("holdNonAdapter").get(0));
		}
		if(info.containsKey("recursiveAdapterTrimming")){
			recursiveAdapterTrimming = Boolean.parseBoolean(info.get("recursiveAdapterTrimming").get(0));
		}
		if(info.containsKey("writeNonAdapter")){
			writeNonAdapter = Boolean.parseBoolean(info.get("writeNonAdapter").get(0));
		}
//		if(info.containsKey("adapterTrimmed")){
//			adapterTrimmed = Boolean.parseBoolean(info.get("adapterTrimmed").get(0));
//		}
		if(info.containsKey("guessAdapter")){
			guessAdapter = Boolean.parseBoolean(info.get("guessAdapter").get(0));
			adapterTrimmed = false;
			if(guessAdapter == false){

			}
			else{
				if( inputType == null && species==null){
					IO.warning("if guessAdapter= is used without genome assembly (species= ), "
							+ "then the type of the input data must be specified inputType=[fastq,fasta,rc]");
					System.exit(1);
				}
				if( species == null){
					IO.writeToCommandLineL1("Will try to guess adapter sequence without genome sequene (NOT RECOMMENDED). Input type is: "+inputType);
					IO.log(logFile, 1, "Will try to guess adapter sequence without genome sequene. Input type is: "+inputType, true);
					String adapterT = Preproc.getMostLikelyAdapter(input, adapterMinLength, 
							17,inputType);
					adapter = adapterT;
					IO.writeToCommandLineL1("Most likely adapter: "+adapterT);
					IO.log(logFile, 2, "Most likely adapter: "+adapterT, true);
				}
			}
		}
		if(info.containsKey("adapter")){
			adapter = info.get("adapter").get(0).replaceAll("U", "T");
			adapterTrimmed = false;
			guessAdapter = false;
		}

		/////////////////////////////////////////////////////////////////////////////////
		// MAPPING PARAMETERS
		
		if(info.containsKey("noMM")){
			mm = Integer.parseInt(info.get("noMM").get(0));
		}
		if(info.containsKey("mm")){
			mm = Integer.parseInt(info.get("mm").get(0));
		}
		if(info.containsKey("seed")){
			seed = Integer.parseInt(info.get("seed").get(0));
			if(seed <= 10){
				IO.warning("Seed parameter is out of acceptable range (>= 10)!! \n"
						+ "       Will quit now!");
				System.exit(0);
				
			}
		}
		if(info.containsKey("alignType")){
			alignType = info.get("alignType").get(0);
		}
		if(info.containsKey("bowtieReportType")){
			bowtieReportType = info.get("bowtieReportType").get(0);
		}
		if(info.containsKey("bowtieReportCount")){
			if(bowtieReportType.equals("-a")){
				IO.warning("bowtieReportCount cannot have any value if bowtieReportType is '-a'. Will ignore this parameter.");
				bowtieReportCount = "";
			}

			else{
				bowtieReportCount = info.get("bowtieReportCount").get(0);
			}
		}
		if(info.containsKey("bowtieAdd")){
			bowtieAdd = info.get("bowtieAdd").get(0);
		}



		if(info.containsKey("libAlignType")){
			libAlignType = info.get("libAlignType").get(0);
		}
		if(info.containsKey("libSeed")){
			libSeed = Integer.parseInt(info.get("libSeed").get(0));
		}

		if(info.containsKey("chunkmbs")){
			chunkmbs = Integer.parseInt(info.get("chunkmbs").get(0));
		}
		if(info.containsKey("p")){
			p = Integer.parseInt(info.get("p").get(0));
		}
		
		if(info.containsKey("direct")){
			direct = Boolean.parseBoolean(info.get("direct").get(0));
		}
		
				
		if(info.containsKey("solid")){
			if(Boolean.parseBoolean(info.get("solid").get(0))){
				IO.writeToCommandLineL1("Input is SOLiD data");
				IO.log(GVars.logFile, 1, "Input is solid", true);
				colorIndex = "_C";
				colorFlag = "-C --col-keepends";
				solid=true;
			}
		}

		if(info.containsKey("microRNAmappingOrientation")){
			microRNAmappingOrientation = info.get("microRNAmappingOrientation").get(0);
			
			if(microRNAmappingOrientation.equals("--nofw") || microRNAmappingOrientation.equals("--norc")){
				
			}
			else{
				IO.warning("microRNAmappingOrientation needs to be either --norc (only forward strand) or --nofw (only reverse strand)!  "
						+ "Will quit now. ");
				System.exit(1);
			}
		}

		if(info.containsKey("libsmappingOrientation")){
			libsmappingOrientation = info.get("libsmappingOrientation").get(0);
			if(libsmappingOrientation.equals("--nofw") || libsmappingOrientation.equals("--norc")){
				
			}
			else{
				IO.warning("libsmappingOrientation needs to be either --norc (only forward strand) or --nofw (only reverse strand)!  "
						+ "Will quit now. ");
				System.exit(1);
			}
		}

		if(info.containsKey("tRNAmappingOrientation")){
			tRNAmappingOrientation = info.get("tRNAmappingOrientation").get(0);
			if(tRNAmappingOrientation.equals("--nofw") || tRNAmappingOrientation.equals("--norc")){
				
			}
			else{
				IO.warning("tRNAmappingOrientation needs to be either --norc (only forward strand) or --nofw (only reverse strand)!  "
						+ "Will quit now. ");
				System.exit(1);
			}
		}


		if(info.containsKey("fileType")){
			fileType = info.get("fileType").get(0);
		}
		

		
		/////////////////////////////////////////////////////////////////////
		/// LIBRARIES
		if(info.containsKey("microRNA")){
			microRNA = info.get("microRNA").get(0);
		}
		if(info.containsKey("homolog")){
			homolog = info.get("homolog").get(0);
		}
		if(info.containsKey("hairpin")){
			hairpin = info.get("hairpin").get(0);
		}
		if(info.containsKey("mature")){
			mature = info.get("mature").get(0);
		}
		if(info.containsKey("rnaSeqLib")){
			rnaSeqLib = info.get("rnaSeqLib").get(0);
		}
		else{
			if(RNAseq){
				IO.warning("rnaSeqLib= needs to be specified for (simple) RNA-seq analysis. Will quit now. ");
				System.exit(1);
			}
		}

		if(info.containsKey("calcFeaturesMiR")){
			calcFeaturesMiR = Boolean.parseBoolean(info.get("calcFeaturesMiR").get(0));
		}

		if(info.containsKey("calcFeaturesDominant")){
			calcFeaturesDominant = Boolean.parseBoolean(info.get("calcFeaturesDominant").get(0));
		}
		
		
		
		///////////////////////////////////////////////////////////////////////////////////////////
		/// preprocessing 
		if(info.containsKey("minReadLength")){
			minReadLength = Integer.parseInt(info.get("minReadLength").get(0));
		}
		if(info.containsKey("maxReadLength")){
			maxReadLength = Integer.parseInt(info.get("maxReadLength").get(0));
		}
		if(info.containsKey("minRC")){
			minRC = Integer.parseInt(info.get("minRC").get(0));
		}
		if(info.containsKey("removeBarcode")){
			removeBarcode = Integer.parseInt(info.get("removeBarcode").get(0));
		}
		if(info.containsKey("remove3pBases")){
			remove3pBases = Integer.parseInt(info.get("remove3pBases").get(0));
		}
		if(info.containsKey("sep")){
			sep = info.get("sep").get(0);
		}

		if(info.containsKey("iterative5pTrimming")){
			iterative5pTrimming = Integer.parseInt(info.get("iterative5pTrimming").get(0));
		}
		
		
		
		//////////////////////////////////////////////////
		///// isoMiR parameters
		if(info.containsKey("isomiRseed")){
			isomiRseed = Integer.parseInt(info.get("isomiRseed").get(0));
		}
		if(info.containsKey("isoMiR")){
			isoMiR = Boolean.parseBoolean(info.get("isoMiR").get(0));
		}
		if(info.containsKey("isoLibs")){
			isoLibs = Boolean.parseBoolean(info.get("isoLibs").get(0));
		}
		if(info.containsKey("minRCiso")){
			minRCiso = Integer.parseInt(info.get("minRCiso").get(0));
		}
		
		if(info.containsKey("fullIsoStat")){
			fullIsoStat = Boolean.parseBoolean(info.get("fullIsoStat").get(0));
		}

		

		
		//////////////////////////////////////////////////////////////////
		///// OUTPUT OPTIONS
		if(info.containsKey("graphics")){
			graphics = Boolean.parseBoolean(info.get("graphics").get(0));
		}
		if(info.containsKey("plotLibs")){
			plotLibs = Boolean.parseBoolean(info.get("plotLibs").get(0));
		}
		if(info.containsKey("plotMiR")){
			plotMiR = Boolean.parseBoolean(info.get("plotMiR").get(0));

		}
		if(info.containsKey("plotSec")){
			plotSec = Boolean.parseBoolean(info.get("plotSec").get(0));
			calcFeaturesMiR = true;

		}
		
		if(info.containsKey("bedGraph")){
			bedGraph = Boolean.parseBoolean(info.get("bedGraph").get(0));
		}
		if(info.containsKey("bedGraphMode")){
			bedGraphMode = info.get("bedGraphMode").get(0);
		}
		if(info.containsKey("bedGraphIntervals")){
			bedGraphIntervals = info.get("bedGraphIntervals").get(0);
		}
		
		
		if(info.containsKey("minRCplotLibs")){
			minRCplotLibs = Integer.parseInt(info.get("minRCplotLibs").get(0));
		}
		if(info.containsKey("minRCplotMiR")){
			minRCplotMiR = Integer.parseInt(info.get("minRCplotMiR").get(0));
		}
		if(info.containsKey("maxLenSecStruc")){
			maxLenSecStruc = Integer.parseInt(info.get("maxLenSsomecStruc").get(0));
		}
		if(info.containsKey("maxLenPlot")){
			maxLenPlot = Integer.parseInt(info.get("maxLenPlot").get(0));
		}

		if(info.containsKey("libsStringTypes")){
			libsStringTypes =  info.get("libsStringTypes").get(0);
			String[] f1 = libsStringTypes.split("\\|");
			if(info.containsKey("libsStringNames")){
				libsStringNames =  info.get("libsStringNames").get(0);
				String[] f2 = libsStringNames.split("\\|");
				if(f1.length != f2.length){
					IO.warning("libsStringType and libsStringNames must contain the same number of elements separated by | ");
					System.exit(1);
				}
					
			}
			else{
				IO.warning("If you set libsStringType you need to set libsStringNames as well! ");
				System.exit(1);
			}
		}

		//////////////////////////////////////////////////////////////////
		///// MAKE GENOME DISTRIBUTION
		
		if(info.containsKey("assocFile")){
			assocFile = info.get("assocFile").get(0);
		}
		if(info.containsKey("splitToSpecies")){
			splitToSpecies = Boolean.parseBoolean(info.get("splitToSpecies").get(0));
		}
		
		if(info.containsKey("mainSpecies")){
			mainSpecies = info.get("mainSpecies").get(0);
		}

		if(info.containsKey("chromosomeLevel")){
			chromosomeLevel = Boolean.parseBoolean(info.get("chromosomeLevel").get(0));
		}

		if(info.containsKey("genomeDistunique")){
			genomeDistunique = Boolean.parseBoolean(info.get("genomeDistunique").get(0));
		}
		
		if(info.containsKey("chrMappingByLength")){
			chrMappingByLength = Boolean.parseBoolean(info.get("chrMappingByLength").get(0));
		}

		if(info.containsKey("chromosomes")){
			chromosomes = info.get("chromosomes").get(0);
		}
		ap.add("chrLength");
		
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//// expression profiling parameters

		if(info.containsKey("winUpMir")){
			winUpMir = Integer.parseInt(info.get("winUpMir").get(0));
		}
		if(info.containsKey("winDownMir")){
			winDownMir = Integer.parseInt(info.get("winDownMir").get(0));
		}
		if(info.containsKey("winUpLibs")){
			winUpLibs = Integer.parseInt(info.get("winUpLibs").get(0));
		}
		if(info.containsKey("winDownLibs")){
			winDownLibs = Integer.parseInt(info.get("winDownLibs").get(0));
		}		
		
		if(info.containsKey("winUpPre")){
			winUpPre = Integer.parseInt(info.get("winUpPre").get(0));
		}
		if(info.containsKey("winDownPre")){
			winDownPre = Integer.parseInt(info.get("winDownPre").get(0));
		}	

		if(info.containsKey("hierarchical")){
			hierarchical = Boolean.parseBoolean(info.get("hierarchical").get(0));
		}
		
		//################
		// novel microRNAs


		
		if(info.containsKey("kingdom")){
			kingdom = info.get("kingdom").get(0);
		}

		if(info.containsKey("predict")){
			predict = Boolean.parseBoolean(info.get("predict").get(0));
			
			if(GVars.species == null){
				IO.warning("When predict=true, then a genome assembly must be given with species=. \n"
						+ "         Will quit now!");
				System.exit(0);
			}
			
			if(kingdom.equals("animal")){
				maxDistNovel = maxDistAnimal;
				novelHomolog = "ssa:hme:ssc:pma:dmo:hhv6b:lmi:hiv1:pmi:der:tur:egr:ccr:dev:cte:dvi:sko:lva:rlcv:ssy:sci:cla:mghv:api:ocu:bbe:ola:"
						+ "hvt:fru:sla:sv40:tni:ipu:str:odi:rrv:dwi:kshv:aqu:aae:sma:sme:cel:bkv:pol:bta:dgr:lgi:dpe:smr:aja:jcv:cfa:ppa:mdv2:"
						+ "mdv1:ppc:tgu:mja:dps:dpu:hhi:bhv5:dya:blv:bhv1:nlo:pxy:bma:ebv:ppy:isc:aca:mse:gga:eca:asu:bmo:bpcv2:bpcv1:mcmv:ggo:cgr:"
						+ "rmi:dan:dre:nve:nvi:mcv:prd:chi:tre:bfl:oha:hru:sha:gpy:rno:bfv:prv:pbi:hsa:dse:ame:mdo:spu:dsi:hbv:xbo:cqu:cin:lca:emu:mml:"
						+ "ngi:tca:mmu:hcmv:iltv:hco:tch:lco:meu:oan:crm:lla:oar:mne:ddi:sja:xla:hvsa:cbn:ptr:csa:cbr:gsa:efu:hsv2:hsv1:xtr:aga:dme:hma:age";
			}
			else if(kingdom.equals("plant")){
				maxDistNovel = maxDistPlant;
				novelHomolog="cca:gso:sof:lus:ccl:ssl:vun:zma:ata:bra:pab:ssp:egu:ath:bna:han:pin:cpa:far:har:tae:hex:atr:pra:ghb:ctr:peu:lja:cln:hvu:vvi:mtr:pvu:"
						+ "aly:ghr:aqc:ahy:ama:stu:mdm:hbr:gma:cme:amg:bol:sly:psj:hci:pgi:smo:aau:bgy:tcc:rgl:gra:cre:hpa:gar:bcy:nta:mes:pta:hpe:ptc:rco:crt:"
						+ "pti:ppe:bdi:dpr:htu:pde:ttu:ppt:csi:esi:sbi:osa";
			}
			else{
				IO.warning("kingdom= must be either plant or animal");
				System.exit(1);
			}
			
		}
		if(info.containsKey("flankNovel")){
			flankNovel = Integer.parseInt(info.get("flankNovel").get(0));
		}
		if(info.containsKey("maxDistNovel")){
			maxDistNovel = Integer.parseInt(info.get("maxDistNovel").get(0));
		}
		if(info.containsKey("novelName")){
			novelName = info.get("novelName").get(0);
		}
		if(info.containsKey("novelHomolog")){
			novelHomolog = info.get("novelHomolog").get(0);
		}
		
		if(info.containsKey("startCountNovel")){
			startCountNovel = Integer.parseInt(info.get("startCountNovel").get(0));
		}
		

		if(info.containsKey("chrMappingByLength")){
			chrMappingByLength = Boolean.parseBoolean(info.get("chrMappingByLength").get(0));
		}

		
		
	}

}
