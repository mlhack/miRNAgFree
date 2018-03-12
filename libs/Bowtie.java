package libs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipOutputStream;

import sequences.SeqUtil;



/**
 * Functions related to read alignment by Bowtie 1
 * 
 *
 */
public class Bowtie {

	
	/**
	 * Convert bowtie output to fasta
	 * @param bowtieParsed
	 * @param readsFile
	 */
	public static void makeReadsFasta(String bowtieParsed, String readsFile){
		
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter(readsFile,true));
			BufferedReader reader = new BufferedReader ( new FileReader(bowtieParsed));
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\\t");
				String seq = f[4];
				if(f[1].equals("-")){
					seq = SeqUtil.getReverseComplementarySequence(f[4]);
				}
				writer.write(">"+f[0]+"\n");
				writer.write(seq+"\n");
			}
			reader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	

	/**
	 * 
	 * @param infile
	 * @param outfile
	 * @param parse: parse the output file (i.e. try to get the best mappings)
	 * @param setParameters: the mapping parameters are set to those given by the species name like hg19#2#19#v  
	 * @param highRedundant: if highly redundant reads are mapped
	 * @param removeMapped: the reads mapped are removed from the infile (!!! when annotations are mapped, this needs to be set to false !!!)
	 * @return
	 */
	public static AlignData genomeAlignBowtie(String infile,String outfile, boolean parse, boolean setParameters, boolean highRedundant, boolean removeMapped ){

		String[] species  = GVars.species.split(":");
		File fn = new File(outfile);
		String name = Util.getFileBaseName(fn.getName());
		String parsed = fn.getParent()+File.separator+name+".parsed";
		if(species.length == 1){

			IO.log(GVars.logFile, 1, "Will use "+species.length+" species for alignment", parse);;
			if(setParameters)
				setMappingParamters(species[0]);
			
			if(highRedundant){
				GVars.tempBowtieReportType = "-k";
				GVars.tempBowtieReportCount = "2";
			}
			
			genomeAlign(infile,GVars.index+File.separator+GVars.species+GVars.colorIndex,
					outfile,parsed,true);
			
			
		}
		else{

			IO.log(GVars.logFile, 1, "Will use "+species.length+" species for alignment", parse);;
			new File(outfile).delete();
			for(String specie : species){

				if(setParameters)
					setMappingParamters(specie);

				if(highRedundant){
					GVars.tempBowtieReportType = "-k";
					GVars.tempBowtieReportCount = "2";
				}
				
				String[] spec = specie.split("#");
				String specieStr = spec[0];

				genomeAlign(infile,GVars.index+File.separator+specieStr+GVars.colorIndex, 
						GVars.output+File.separator+specie+".txt",parsed,false);	

				IO.copy(GVars.output+File.separator+specie+".txt", 
						outfile, true);
				if(GVars.remove)
					new File(GVars.output+File.separator+specie+".txt").delete();

			}
			// sort by reads
			SortText.sortAsc(outfile, 0, false);
		
			
		if(GVars.alignType.equals("n"))
			Bowtie.parseGenomeAlignment(outfile, parsed,GVars.tempSeed);
		else
			Bowtie.parseGenomeAlignment(outfile,parsed);


		
		}

		if(!(new File(parsed).isFile()))
			return null;
	
		if(removeMapped){
			Set<String> mappedID = getBowtieReferences(parsed, 0);
			Util.removeReadsFasta(infile,mappedID);
			Util.sortSRNAbenchFormat(infile);
		}
		AlignData data = Stat.profileFromBowtieOut(parsed);
		data.library = name;
	

		return data;
	}
	
	/**
	 * gives back a set with all reference names that are in a bowtie output file
	 * @param bowtieFile
	 * @index: the column name: 0 --> ids; 2 --> reference names
	 * @return
	 */
	public static Set<String> getBowtieReferences(String bowtieFile,int index){
		Set<String> back = new HashSet<String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(bowtieFile));
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				back.add(f[index]);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return back;
	}
	/**
	 * gives back a set with all reference names that are in a bowtie output file
	 * @param bowtieFile
	 * @index: the column name: 0 --> ids; 2 --> reference names
	 * @return
	 */
	public static Set<String> getMappedReads(String bowtieFile){
		Set<String> back = new HashSet<String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(bowtieFile));
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				String readS = new String(f[4]);
				if(f[0].equals("-"))
					readS = SeqUtil.getReverseComplementarySequence(f[4]);
				
				back.add(readS);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return back;
	}
	
	/**
	 * 
	 * this function launches the Bowtie program for grouped input!!!
	 * 
	 * @param input the input in fasta
	 * @param index the index base file name
	 * @param output
	 * @param parse --> parse the output or not
	 */
	public static void genomeAlign(String input, String index,String output, String parsedFile, boolean parse){

		if(GVars.tempAlingType.equals("n")){

			String parameters = "--chunkmbs "+GVars.chunkmbs+" -p "+GVars.p+" -n "+GVars.tempMM+" -l "+	GVars.tempSeed+" "+GVars.tempBowtieReportType+" "+GVars.tempBowtieReportCount+"  "+GVars.bowtieAdd+"  -e 20000 "+GVars.colorFlag+" "+GVars.fileType+" "+index+" "+input+" "+output;
			IO.log(GVars.logFile, 5, parameters, true);
			if(GVars.direct)
				Exec.bowtieAlignDirect(GVars.bowtie, parameters);
			else
				Exec.bowtieAlign(GVars.bowtie, parameters);
			
			if(parse)
				parseGenomeAlignment(output, parsedFile, GVars.seed);
		}
		else if(GVars.tempAlingType.equals("v")){

			String parameters = "--chunkmbs "+GVars.chunkmbs+" -p "+GVars.p+" -v "+GVars.tempMM+"  "+GVars.tempBowtieReportType+" "+GVars.tempBowtieReportCount+"  "+GVars.bowtieAdd+"  -e 20000  "+GVars.colorFlag+" "+GVars.fileType+" "+index+"   "+input+" "+output;
			IO.log(GVars.logFile, 5, parameters, true);
			if(GVars.direct)
				Exec.bowtieAlignDirect(GVars.bowtie, parameters);
			else
				Exec.bowtieAlign(GVars.bowtie, parameters);
			if(parse)
				parseGenomeAlignment(output, parsedFile);
		}
		else{
			IO.warning("the alignment type needs to be either 'n' or 'v' !!!!!!");
		}
	}


	
	/**
	 * parses and original bowtie output file --> holds only the best reads seed extension & number of mismatches outside the seed
	 * @param file
	 * @param parsedFile: the output file where the parsed file will be written
	 * @param seedLength
	 * @return the filename of the parsed file. Note: the string can be null
	 */
	public static boolean parseGenomeAlignment(String file, String parsedFile, int seedLength ){


		File f1 = new File(file);
		if(!f1.exists()){
			IO.log(GVars.logFile, 4, "No mapping file has been found", true);
			return false;
		}
		
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter(parsedFile));
			BufferedReader reader = new BufferedReader (new FileReader(file));
			String line = reader.readLine();
			String oldID = null;
			List<String> list = new ArrayList<String>();
			if(line != null){
				oldID = line.split("\t")[0];
				list.add(line);
			}
			while((line = reader.readLine()) != null){

				String[] f = line.split("\t");
				if(f[0].equals(oldID)){
					list.add(line);
				}
				else{
					list = getBestSeedAlignments(list, seedLength);
					List<String> best = getBestAlignments(list, seedLength);
					writeBowtieParsed(writer, best);						
					// add
					oldID=f[0];
					list = new ArrayList<String>();
					list.add(line);
				}
				
			}
			// write out last block
			list = getBestSeedAlignments(list, seedLength);
			List<String> best = getBestAlignments(list, seedLength);
			if(best.size() == 0){
				IO.warning("empty set of reads, maybe due to bad format!");
			}
			else
				writeBowtieParsed(writer, best);						

			writer.close();
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return true;
	}

	
	/**
	 * parses and original bowtie output file --> holds only the best reads seed extension & number of mismatches outside the seed
	 * writes and output file with the 'file' name + 'parsed' extension
	 * @param file
	 * @return
	 */
	public static boolean parseGenomeAlignment(String file, String parsedFile){

		
		File f1 = new File(file);
		if(!f1.exists()){
			IO.log(GVars.logFile, 4, "No mapping file has been found", true);
			return false;
		}
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter(parsedFile));

			BufferedReader reader = new BufferedReader (new FileReader(file));
			String line = reader.readLine();
			String oldID = null;
			List<String> list = new ArrayList<String>();
			int sl = 0;
			if(line != null){
				oldID = line.split("\t")[0];
				list.add(line);
				sl = line.split("\t")[4].length();
			}
			while((line = reader.readLine()) != null){

				String[] f = line.split("\t");
				if(f[0].equals(oldID)){
					list.add(line);
				}
				else{
					list = getBestSeedAlignments(list, sl);
					List<String> best = getBestAlignments(list, sl);
					writeBowtieParsed(writer, best);						
					// add
					oldID=f[0];
					sl = f[4].length();
					list = new ArrayList<String>();
					list.add(line);
				}
				
			}
			// write out last block
			list = getBestSeedAlignments(list, sl);
			List<String> best = getBestAlignments(list, sl);
			writeBowtieParsed(writer, best);						

			writer.close();
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return true;
	}
	
	/**
	 * write out a manipulated bowtie output file
	 * @param writer
	 * @param list
	 * @param count
	 */
	public static void writeBowtieParsed(BufferedWriter writer,List<String> list){
		
		for(String line : list){
			String[] f = line.split("\t");
			try {
				writer.write(f[1]+"\t"+f[2]+"\t"+f[3]+"\t"+f[4]+"\t"+f[5]+"\t"+f[6]+"\t"+list.size()+"\t"+f[0]);
				if(f.length == 8){
					writer.write("\tNA");
				}
				else{
					writer.write("\t"+f[8]);
				}
				writer.write("\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	/**
	 * gives back a list with the best seed alignments --> this is important when bowtie output files from different runs are combined
	 * @param bowtieList
	 * @param seedLength
	 * @return
	 */
	public static List<String> getBestSeedAlignments(List<String> bowtieList,int seedLength){
		
		List<String> back = new ArrayList<String>();
		int bestSeed = 100000;
		for(String line : bowtieList){

			String[] f = line.split("\t");
			// no mismatches --> add line directly		
			if(f.length == 7){
				if(bestSeed > 0)
					back = new ArrayList<String>();
				
				bestSeed = 0;
				back.add(line);
			}
			else if(f.length < 7){
				IO.warning("Error in input line: "+line+". This should be a bowtie1 output line!");
			}
			else{
				boolean checkMM = checkMismatchString(f[7]);
				if(checkMM){
					int seedMM = getSeedMM(f[7], seedLength);
					if(seedMM == bestSeed){
						back.add(line);
					}
					else if(seedMM < bestSeed){
						back = new ArrayList<String>();
						back.add(line);
						bestSeed = seedMM;
					}
				}
			}

		}
		return back;

	}
	
	/**
	 * filter the best alignments --> do a seed extension
	 * @param bowtieList
	 * @param seed
	 * @return
	 */
	private static List<String> getBestAlignments(List<String> bowtieList,int seedLength){
		
		List<String> back = new ArrayList<String>();
		int bestLen = -1;
		int bestNoSeed = 1000000;

		for(String line : bowtieList){

			String[] f = line.split("\t");
			// no mismatches --> add line directly		
			if(f.length == 7){
				bestLen = f[4].length();
				bestNoSeed = 0;
				back.add(0+"\t"+line);
			}
			else{
				boolean checkMM = checkMismatchString(f[7]);

				if(checkMM){
					int tempLen = getBestLength(f[7], seedLength, f[4].length());
					int nonSeedMM = getNonSeedMM(f[7], seedLength);
					if(tempLen == bestLen && nonSeedMM == bestNoSeed){
						back.add(nonSeedMM+"\t"+line);
					}
					else if(tempLen > bestLen && nonSeedMM == bestNoSeed){
						back = new ArrayList<String>();
						back.add(nonSeedMM+"\t"+line);
						bestLen = tempLen;
					}
					else if(tempLen == bestLen && nonSeedMM < bestNoSeed){
						back = new ArrayList<String>();
						back.add(nonSeedMM+"\t"+line);
						bestNoSeed = nonSeedMM;					
					}
					else if(tempLen > bestLen && nonSeedMM < bestNoSeed){
						back = new ArrayList<String>();
						back.add(nonSeedMM+"\t"+line);
						bestNoSeed = nonSeedMM;					
						bestLen = tempLen;
					}
				}
				else{
					IO.warning("BAD LINE: "+line);
				}
			}

		}
		return back;
	}
	
	static boolean checkMismatchString(String misMatchString){
		String[] code = misMatchString.split(",");
		
		for(String c : code){
			String t = c.split(":")[0];
			try{
			int l = Integer.parseInt(t);
			}catch(NumberFormatException e){
				return false;
			}
		}
		return true;
	}

	/**
	 * gets back the longest alignment surpassing the bowtie alignment seed
	 * @param misMatchString
	 * @param seedLength
	 * @param readLength the length of the read --> if no mismatch out of the seed is found --> give back the read length
	 * @return
	 */
	static int getBestLength(String misMatchString, int seedLength, int readLength){
		

		String[] code = misMatchString.split(",");
		
		for(String c : code){
			String t = c.split(":")[0];
			int l = Integer.parseInt(t);
			if(l + 1 > seedLength ){
				return l ;
			}				
		}
		return readLength;
	}
	
	/**
	 * detects the number of mismatches within the seed region
	 * @param misMatchString
	 * @param seedLength
	 * @return
	 */
	static int getSeedMM(String misMatchString, int seedLength){

		String[] code = misMatchString.split(",");
		int count = 0;
		for(String c : code){
			String t = c.split(":")[0];
			try{
				int l = Integer.parseInt(t);
				if(l + 1 <= seedLength ){
					count++;
				}			
			}
			catch(NumberFormatException e){
				IO.warning("Number format exception for mismatch string: "+misMatchString+" in getSeedMM function");
			}
		}
		return count;

	}
	
	
	/**
	 * detects the number of mismatches outside the seed regio
	 * @param misMatchString
	 * @param seedLength
	 * @return
	 */
	static int getNonSeedMM(String misMatchString, int seedLength){

		String[] code = misMatchString.split(",");
		int count = 0;
		for(String c : code){
			String t = c.split(":")[0];
			int l = Integer.parseInt(t);
			if(l + 1 > seedLength ){
				count++; ;
			}				
		}
		return count;

	}
	
	/**
	 * this function sets the temporary mapping parameters given as extension to the library name. !! it does not reset the temporary mapping parameters !!!
	 * @param libName
	 */
	public static void setMappingParamters(String libName){

		
		////////////////////////////////////////////////
		/// parse out the parameters given in the string 
		String[] par = libName.split("#");
		if(par.length > 1){
			GVars.tempMM = Integer.parseInt(par[1]);
			if(par.length > 2){
				GVars.tempSeed = Integer.parseInt(par[2]);
				if(par.length > 3){
					GVars.tempAlingType = par[3];
					if(par.length > 4){
						GVars.tempBowtieReportType = par[4];
						if(par[4].equals("-a")){
							GVars.tempBowtieReportCount = "";
						}
						if(par.length > 5 && !(par[4].equals("-a"))){
							GVars.tempBowtieReportCount = par[5];
						}
					}

				}
			}
		}
		////////////////////////////////////////////////////
	}
	
	/**
	 * Reset the bowtie mapping parameters to default
	 */
	public static void resetMappingParameters(){
		
		GVars.tempAlingType = GVars.alignType;
		GVars.tempBowtieReportCount = GVars.bowtieReportCount;
		GVars.tempBowtieReportType = GVars.bowtieReportType;
		GVars.tempMM = GVars.mm;
		GVars.tempSeed = GVars.seed;
		GVars.mappingOrientation = "";
	}
	
	
	/**
	 * Reads a bowtie output alignment file and gets the most likely adapter back
	 * @param bowtieFile
	 * @param minAdapterLen: the minimum adapter length
	 * @param seedLength: the seed length (to know where the putative adapter starts)
	 * @return
	 */
	public static String getAdapterSequenceFromAlignment(String bowtieFile, int minAdapterLen, int seedLength){
		
		Map<String,int[]> countMap = new Hashtable<String,int[]>();
		try {
			BufferedReader reader = new BufferedReader ( new FileReader(bowtieFile));
			
			String line = null;
			while((line = reader.readLine()) != null){
				
				String[] f = line.split("\t");
				String seq = f[4];
				if(f[1].equals("-"))
					seq = SeqUtil.getReverseComplementarySequence(f[4]);
				
				String adapterFrag =  getPutativeAdapter(f[8], seq, seedLength, minAdapterLen);
				if(adapterFrag == null)
					continue;
				
				if(countMap.containsKey(adapterFrag)){
					countMap.get(adapterFrag)[0]++;
				}
				else{
					int[] t = new int[1];
					t[0]++;
					countMap.put(adapterFrag, t);
				}
			}
			
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String back = "XXXXX";
		int highest = 0;
		for(String key : countMap.keySet()){
			if(SeqUtil.ratioOfMostFrequentBase(key) < 0.7 && !(key.contains("N"))){
				if(countMap.get(key)[0] > highest){
					back = key;
					highest = countMap.get(key)[0];
				}
			}
		}
		return new String(back);
		
	}



	/**
	 * 
	 * @param outPath
	 * @param libName
	 * @param map
	 * @param totalReadsLib
	 * @param totalReads
	 */
	public static String[] writeBowtieMapOut(String outPath, String libName, Map<String,double[]> map, int rcSense, int rcASense, int totalReads){

		String senseOut = outPath + File.separator + libName+"_sense.grouped";
		String aSenseOut = outPath + File.separator + libName+"_antisense.grouped";
		try {


			BufferedWriter writerS = new BufferedWriter( new FileWriter (senseOut));
			BufferedWriter writerAS = new BufferedWriter( new FileWriter (aSenseOut));

			writerS.write(Util.getGroupedHeader()+"\n");
			writerAS.write(Util.getGroupedHeader()+"\n");
			
			for(String name : map.keySet()){
				
				// this is the count for sense
				if(map.get(name)[0] > 0){
					double rpmlib = 1000000*map.get(name)[1]/(double)rcSense;
					double rpmall = 1000000*map.get(name)[1]/(double)totalReads;
					writerS.write(name+"\t"+map.get(name)[0]+"\t"+map.get(name)[1]+"\t"+map.get(name)[2]+"\t"+rpmlib+"\t"+rpmall+"\n");
				}
				// count for an
				if(map.get(name)[3] > 0){
					double rpmlib = 1000000*map.get(name)[4]/(double)rcASense;
					double rpmall = 1000000*map.get(name)[4]/(double)totalReads;
					
					writerAS.write(name+"\t"+map.get(name)[3]+"\t"+map.get(name)[4]+"\t"+map.get(name)[5]+"\t"+rpmlib+"\t"+rpmall+"\n");
				}

			}
			writerS.close();
			writerAS.close();
		}catch (FileNotFoundException e){
			IO.log(GVars.logFile, 4, "File not found exception in Bowtie.writeBowtieMapOut ", true);
			return null;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Error while writing to file in Bowtie.writeBowtieMapOut ", true);
			return null;
		}
		Sort.sortListBigger2Smaller(senseOut, 2, true);
		Sort.sortListBigger2Smaller(aSenseOut, 2, true);
		String[] back = new String[2];
		back[0] = senseOut;
		back[1] = aSenseOut;
		return back;


	}

	/**
	 * 
	 * @param outPath
	 * @param libName
	 * @param map
	 * @param totalReadsLib
	 * @param totalReads
	 */
	public static void writeBowtieMapOutRNAseq(String outfile,  Map<String,double[]> map, AlignData ad, 
			Map<String,Integer> referenceLengthMap, String strand){

		try {


			BufferedWriter writerS = new BufferedWriter( new FileWriter (outfile));


			writerS.write(Util.getGroupedHeader()+"\n");


			for(String name : map.keySet()){

				// this is the count for sense
				if(strand.equals("sense")){
					if(map.get(name)[0] > 0){
						double rpmlib = 1000000d*map.get(name)[1]/((double)ad.rcS*(double)referenceLengthMap.get(name));
						writerS.write(name+"\t"+map.get(name)[0]+"\t"+map.get(name)[1]+"\t"+map.get(name)[2]+"\t"+rpmlib+"\t"+0+"\n");
					}
				}
				else{
					if(map.get(name)[3] > 0){
						double rpmlib = 1000000d*map.get(name)[4]/((double)ad.rcAS*(double)referenceLengthMap.get(name));
						writerS.write(name+"\t"+map.get(name)[3]+"\t"+map.get(name)[4]+"\t"+map.get(name)[5]+"\t"+rpmlib+"\t"+0+"\n");
					}	
				}


			}
			writerS.close();

		}catch (FileNotFoundException e){
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	/**
	 * This function: i) mapps the input reads against an index ('mappingOrientation' parameter!), ii) determines the number of mapped reads 
	 * @param outDir:
	 * @param index
	 * @param name: the name of the library
	 * @param input
	 * @param colorFlag
	 * @param colorIndex
	 * @param totalReads
	 * @return: the name of the parsed file. Can return 'null'
	 */
	public static String mapToIndex(String outDir, String index, String name, String input, String colorFlag, String colorIndex, boolean parse){
		
		String output = outDir+File.separator+name+".bowtieOut";
		String parsed = outDir+File.separator+name+".parsed";
		String parameters = "";
		
		
		boolean check = true;
		if(GVars.tempAlingType.equals("v")){
			parameters = "--chunkmbs "+GVars.chunkmbs+" "+GVars.fileType+" -p "+GVars.p+" -v "+GVars.tempMM+"  "
					+ ""+GVars.tempBowtieReportType+" "+GVars.tempBowtieReportCount+" "+GVars.mappingOrientation+"   "
							+ ""+GVars.bowtieAdd+" -e 200000 "+colorFlag+" " + index+colorIndex+"  "+input+" "+output;
			IO.log(GVars.logFile, 6, parameters, true);

//			Exec.bowtieAlign(GVars.bowtie, parameters);
			
			if(GVars.direct)
				Exec.bowtieAlignDirect(GVars.bowtie, parameters);
			else
				Exec.bowtieAlign(GVars.bowtie, parameters);
			
			
			if(parse)
				check = Bowtie.parseGenomeAlignment(output, parsed);	
		}
		else if (GVars.tempAlingType.equals("n")){
			parameters = "--chunkmbs "+GVars.chunkmbs+" "+GVars.fileType+" -p "+GVars.p+" -n "+GVars.tempMM+" -l "+GVars.tempSeed+" "+GVars.tempBowtieReportType+" "+GVars.tempBowtieReportCount+
			" "+GVars.mappingOrientation+" "+GVars.bowtieAdd+" -e 200000 " +colorFlag+" " + index+colorIndex+"  "+input+" "+output;
			IO.log(GVars.logFile, 6, parameters, true);

//			Exec.bowtieAlign(GVars.bowtie, parameters);
			
			if(GVars.direct)
				Exec.bowtieAlignDirect(GVars.bowtie, parameters);
			else
				Exec.bowtieAlign(GVars.bowtie, parameters);
			
			if(parse)
				check = Bowtie.parseGenomeAlignment(output, parsed, GVars.tempSeed);	
			
		}
		if(check && parse){
//			new File(output).delete();
			return parsed;
		}
		else if(!(parse)){
			return output;
		}
		else{
			return null;
		}

	}

	/**
	 * 
	 * @param library: the library input string (sRNAbench format)
	 * @param readsFile: the reads file in fasta format (sRNAbench format)
	 * @param overwrite: if the library is fasta format --> true: overwrites the bowtie-index, false: does not overwrite the bowtie-index
	 * @param totalReads
	 */
	public static AlignData profileLibrary(String library, String readsFile, String outPath, 
			boolean overwrite, int totalReads, boolean parse, String mappingOrientation){
		
		String[] index = getIndexAndName(library, overwrite);
		if(index == null){
			return null;
		}

		resetMappingParameters();
		setMappingParamters(library);
//		IO.warning("SET parameters for library in profileLibrary");
		GVars.mappingOrientation= mappingOrientation;
		String parsedFile = mapToIndex(outPath,index[0], index[1], readsFile, GVars.colorFlag, GVars.colorIndex, parse);
		resetMappingParameters();
		if(parsedFile != null){
			AlignData ad = Stat.profileFromBowtieOut(parsedFile);
			ad.library = index[1];
			ad.libraryPathBowtie = index[0];
			ad.libraryPathFasta = index[2];
			Map<String,double[]> countMap = Stat.getCountsBowtieOut(parsedFile);
			String[] outFiles = Bowtie.writeBowtieMapOut(GVars.output, index[1], countMap, ad.rcS, ad.rcAS, totalReads);
			
			if(outFiles == null)
				return null;
			
			ad.senseOut = outFiles[0];
			ad.aSenseOut = outFiles[1];
			return ad;
		}
		return null;
	}
	
	
	/**
	 * 
	 * @param file: the library name as given by the user: either only a name (supposed to be in the database) or with full path
	 * @return: [0]: the mode, either annot (the lib is in fasta, bed or gff format) or index (the library is in bowtie index format)
	 */
	public static String[] getIndexAndNameGenome(String file){
		
		String[] back = new String[3];
		File libF = new File(file);

		String mode = "annot";
		// if no extension is given --> map to a given library
		String f[] = libF.getName().split("#")[0].split("\\.");
		String folder = f[0];

		String libraryFile = GVars.libsPath+File.separator+libF.getName().split("#")[0];

		// if libF.getParent == null --> libraries are expected to be in the database
		if(libF.getParent() != null)
			libraryFile = file.split("#")[0];
		
		if(f.length == 1)
			mode = "index";
		
		back[0] = mode;
		back[1] = folder;
		back[2] = new String(libraryFile);

		return back;

	}
	
	/**
	 * This function takes the name of a library and tests: <br>
	 * 	i)  if the name corresponds to a full path (i.e. it is not in the sRNAbenchDB)<br>
	 *  ii) if the library is in fasta format --> it generates an bowtie-index
	 * @param file
	 * @param overwrite: if the bowtie index should be overwritten in case it exists already
	 * @return [0] --> the path to the bowtie index, [1] the name of the library (removing all formating stuff)
	 */
	public static String[] getIndexAndName(String file, boolean overwrite){
		
		String[] back = new String[3];
		File libF = new File(file);

		// if no extension is given --> map to a given library
		String name = libF.getName().split("#")[0];
		String f[] = libF.getName().split("#")[0].split("\\.");
		String folder = f[0];


		
		String libraryFile = GVars.libsPath+File.separator+libF.getName().split("#")[0];

		// if libF.getParent == null --> libraries are expected to be in the database
		if(libF.getParent() != null)
			libraryFile = file.split("#")[0];


		String index = "";
		/**
		 * the library was given as bowtie index
		 */
		if(f.length == 1){
			index = libraryFile;
			
			back[0] = index;
			back[1] = folder;
			back[2] = libraryFile+".fa";
		}
		else{	
			if(!(name.endsWith(".fa") || name.endsWith(".fasta") || name.endsWith(".mfa")  || name.endsWith(".fas")) ){
				IO.warning("The following library "+name+" seems not to be in fasta format (*.fa, *.fasta, *.mfa). Will skip this library.");
				IO.log(GVars.logFile, 3, "The following library "+name+" seems not to be in fasta format (*.fa, *.fasta, *.mfa). Will skip this library.", true);
				return null;
			}
			
			boolean check = Exec.makeBowtieIndex(GVars.bowtieBuild, libraryFile, new File(libraryFile).getParent() + File.separator + folder,overwrite);
			index =  new File(libraryFile).getParent() + File.separator + folder;
			
			back[0] = index;
			back[1] = folder;
			back[2] = libraryFile;
		}
		

		return back;
	}

	/**
	 * 
	 * @param misMatchString
	 * @param read
	 * @param seedLength
	 * @param minAdapterLen
	 * @return
	 */
	static String getPutativeAdapter(String misMatchString, String read,int seedLength, int minAdapterLen){

		if(misMatchString.equals("NA"))
			return null;
		
		String[] code = misMatchString.split(",");
		int pos = -1;
		for(String c : code){
			String t = c.split(":")[0];
			int l = Integer.parseInt(t);
			if(l + 1 > seedLength ){
				pos = l; // the start position of the adapter in 0-based coordinates
				break;
			}				
		}
		// check if whole adapter fragment is in read
		if(pos >= 0 && read.length() - (pos + 1)>= minAdapterLen ){
			return read.substring(pos, pos+minAdapterLen);
		}
		return null;

	}
	


	/**
	 * takes a bowtie output file and splits it into the different reference sequences (chromosomes)
	 * after splitting, the files are packed into a zip file and sorted by chromosome start
	 * @param file: the bowtie output file (*.parsed)
	 * @param zipFile: the name of the output zip file
	 * @return
	 */
	public static Set<String> splitBowtieOutToChromosSorted(String file,String zipFile) {

		Set<String> back = new HashSet<String>();
		if(!(new File(file).exists())){
			return back;
		}
		try {
			
			ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipFile));

			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			line = reader.readLine();
			if(line == null){
				out.close();
				reader.close();
				return back;
			}
			String[] f = line.split("\\t");
			String oldChr = f[2];
			List<String> outList = new ArrayList<String>();
			outList.add(line);
			back.add(oldChr);
			while((line = reader.readLine()) != null){
				f = line.split("\t");
				if(oldChr.equals(f[2])){
					outList.add(line);
				}
				// they are different --> write out
				else{
					Zip.writeToZip(outList, out, oldChr,3);
					outList = new ArrayList<String>();
					oldChr = f[2];
					outList.add(line);
				}
			}

			Zip.writeToZip(outList, out, oldChr,3);
			
			reader.close();
			out.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return back;

	}
	
	public static void genomeGuess(String input,  String index,  int seedLength, int p,String output){
		
		File fn = new File(output);
		String parsed = fn.getParent()+File.separator+fn.getName().split("\\.")[0]+".parsed";
		String parameters = "-p "+p+" -n 0 -l "+	seedLength+"  -k 2  --chunkmbs 512  --best --strata -e 10000 "+GVars.colorFlag+" "+index+"  "+input+" "+output;

//		Exec.bowtieAlign(GVars.bowtie, parameters);
		
		if(GVars.direct)
			Exec.bowtieAlignDirect(GVars.bowtie, parameters);
		else
			Exec.bowtieAlign(GVars.bowtie, parameters);
		
		
		parseGenomeAlignment(output, parsed,seedLength);
	}
	
	

	/**
	 * 
	 * @param libFile
	 * @param outfile
	 * @param parse
	 * @return: a Map chromosome <--> list of BedDataAnnot (annotable BedData class - annotable with MapData objects, i.e. assigning the reads that map to this element)
	 */
	public static Map<String,List<BedDataAnnot>> getBedObjectMapFa(String libFile, String outfile,boolean parse){
		
		File f = new File(libFile);
		String name = Util.getFileBaseName(f.getName());
		
		if(libFile.endsWith("fa") || libFile.endsWith("fasta")){
			IO.log(GVars.logFile, 1,  "Will map the annotations to the genome "+name, true); 
			IO.log(GVars.logFile, 3,  "The results will be FALSE if these sequences are spliced!! "+name, true); 
		
			
			String parsedFile = GVars.output +File.separator+Util.getFileBaseName(outfile)+".parsed";
			
			GVars.tempSeed = GVars.libSeed;
			GVars.tempAlingType = GVars.libAlignType;
			GVars.tempMM = GVars.libsMM;
			GVars.tempBowtieReportType = "-a";
			GVars.tempBowtieReportCount = "";
			genomeAlignBowtie(libFile, outfile, parse, false, false, false);
			resetMappingParameters();

			Map<String,List<BedDataAnnot>> map = BedDataAnnot.bowtie2BED(parsedFile);
			BedDataAnnot.sortStartAsc(map);
//			BedDataAnnot.removeDuplicates(map);
			return map;
		}
		else if(libFile.endsWith("bed")){
			Map<String,List<BedDataAnnot>> back = BedDataAnnot.read(libFile, false);
			return back;
		}
		// get from gtf file
		else if(libFile.contains(".gtf") || libFile.contains(".gff")){
				
				String[] f1 = libFile.split("#");
				String feature = GVars.gtfFeature;
				if(f1.length >= 2){
					feature = f1[1];		
				}
//				System.out.println(f[0]);
				Map<String,List<BedDataAnnot>> back = GtfProcesser.readFile(f1[0], feature, GVars.gtfAtribute);

				return back;
			}
		
		return null;
	}
	
}
