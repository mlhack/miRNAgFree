package libs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;







//import sRNAbench.Vars;
import sequences.SeqUtil;




/**
 * Preprocessing functions for sRNAbench - miRNAgFree - sRNAblast
 * @author Michael
 *
 */
public class Preproc {

	
	static public int readsRaw = 0; // number of raw reads in the analysis
	static public int readsAdapterFound = 0; // number of reads for which the adapter was found
	static public int readsAdapterNotFound = 0; // number of reads for which the adapter was not found
	static public int readsLengthFilteredMin = 0; // number of reads filtered out because to short
	static public int readsLengthFilteredMax = 0; // number of reads filtered out because to long
	static public int reads = 0; // number of reads in the analysis
	static public int readsUnique = 0; // number of unique reads in analysis
	static public int readsQRCfiltered = 0; // number of quality filtered reads
	static public int readsMaxLengthInput = 0; // maximum read length in input
	static public int readsMaxLengthAnalysis = 0; // maximum read length in input
	static public int filteredRC = 0;
	static public int filteredUR = 0;
	static public boolean preprocFinished = false;
	
	/**
	 * This function takes some preprocessing steps
	 * 1) if the input is solid
	 * 2) if the input is SRA
	 * 3) it runs the input() function to i) perform adapter trimming, ii) convert input to internal sRNAbench file
	 */
	public static boolean preprocessing(){
		
		if(!(GVars.doAlignment)){
			return false;
		}
		
		IO.writeToCommandLineBlockOpen("START WITH THE PRE-PROCESSING OF THE READS");
		IO.log(GVars.logFile, 1, "Start with preprocessing of the reads", true);
		////////////////////////////////////////////////////////////////////////////////////////////
		/// preparing the data


		
		/*
		 * Convert SRA to fastq
		 */
		if((GVars.input.endsWith("sra") || (GVars.input.startsWith("SRR")) || GVars.input.startsWith("ERR"))  ){
			
			String[] files = GVars.input.split(":");
			StringBuilder sb = new StringBuilder();
			String sraTemp = null;
			for(String file : files){ //**************************************
					
				File ip = new File(file);
			
				if(ip.getName().endsWith("fastq.gz") || ip.getName().endsWith("fastq")){
					break;
				}
				
				else if(ip.isFile()){
					IO.log(GVars.logFile, 1, "Local SRA input file! ", true);
					IO.writeToCommandLineL1("Local SRA input file! ");
					sraTemp = new String(ip.getAbsolutePath());
				}
				else if(ip.getName().startsWith("SRR")){
					
					IO.log(GVars.logFile, 1, "Remote SRA input file. Will try to download. ", true);
					IO.writeToCommandLineL1("Remote SRA input file. Will try to download.");
					String urlString = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"+ip.getName().substring(0, 6)+"/"+file+"/"+file+".sra";

					IO.saveUrl(GVars.output+File.separator+file+".sra", urlString);
					sraTemp = GVars.output+File.separator+file+".sra";
				}
				
				else if( ip.getName().startsWith("ERR")){
					
					IO.log(GVars.logFile, 1, "Remote SRA input file. Will try to download. ", true);
					String urlString = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/"+ip.getName().substring(0, 6)+"/"+file+"/"+file+".sra";
					IO.writeToCommandLineL1("Remote SRA input file. Will try to download: "+urlString);

					IO.saveUrl(GVars.output+File.separator+file+".sra", urlString);
					sraTemp = GVars.output+File.separator+file+".sra";
				}
				
				else{
					IO.log(GVars.logFile, 4, "Input file has sra extension, but: i) the file was not found locally or ii) it was not a valid SRR identifier. Will quit now. ", true);
					IO.warning("Input file has sra extension, but: i) the file was not found locally or ii) it was not a valid SRR identifier. Will quit now. ");
					System.exit(1);
				}
				IO.log(GVars.logFile, 1, "Converting SRA file to fastq "+ip.getName(), true);
				IO.writeToCommandLineL1("Converting SRA file to fastq "+ip.getName());
				
				Exec.convertSRA(sraTemp, " -O "+GVars.output, "fastq-dump");

				sb.append(GVars.output+File.separator+ip.getName().split("\\.")[0]+".fastq:");

				
			} //********************************
			GVars.input = sb.toString().substring(0, sb.length()-1);
			GVars.origInput = sb.toString().substring(0, sb.length()-1);
			IO.log(GVars.logFile, 1, "After SRA conversation, new input file: "+GVars.input, true);
			IO.writeToCommandLineL1("After SRA conversation, new input file: "+GVars.input);
			
		}

		if(!(new File(GVars.input).isFile()) && !(GVars.input.contains(":")) ){
			IO.log(GVars.logFile, 1, "Remote input file. Will try to download. ", true);
			IO.writeToCommandLineL1("Remote input file. Will try to download.");
			String fn = new File(GVars.input).getName();
			boolean check = IO.saveUrl(GVars.output+File.separator+fn, GVars.input);
			if(!(check)){
			}
			else{
				GVars.origInput = GVars.input;
				GVars.input = GVars.output+File.separator+fn;

			}
		}

		if(GVars.solid  && GVars.guessAdapter){

			IO.warning("guessAdapter does currently not work with SOLiD! Will quit now!");
			IO.log(GVars.logFile, 4, "uessAdapter= does currently not work with SOLiD!", true);
			System.exit(1);

		}
		/**
		 * ONLY IF SOLID
		 */
		if(GVars.solid){
			prepareSolid();
		}

		/*
		 * 
		 */
		else if(GVars.guessAdapter && GVars.species != null){

			IO.writeToCommandLineL1("Start with guess adapter sequence");
			IO.log(GVars.logFile, 1, "Start with guess adapter sequence", true);

			String out = GVars.output+File.separator+"gtemp.fastq";
			//			GVars.writeToCommandLineL1(out);
			IO.makeSubFastq(GVars.input, out, 1000000);

			Bowtie.genomeGuess(out, GVars.index+File.separator+GVars.species.split(":")[0]+GVars.colorIndex, 
					GVars.seed,  GVars.p, GVars.output+File.separator+"genome.txt");
			

			String adapterSeq = Bowtie.getAdapterSequenceFromAlignment(GVars.output+File.separator+"genome.parsed", GVars.adapterMinLength, GVars.seed);
			GVars.adapter = adapterSeq;
			IO.warning("infered adapter sequence: "+adapterSeq);
			IO.log(GVars.logFile, 2, "Infered adapter sequence: "+adapterSeq, true); 
			new File(out).delete();

		}
		///// END PREPARING THE DATA
		/////////////////////////////////////////////////////////////////////////////////////
		
		
		// read the input
		input();
//		IO.warning("input file is now: "+GVars.input);
		/**
		 * i) 
		 */
		if(GVars.iterative5pTrimming > 0){
			iterative5pTrimming();
		}
		
		preprocFinished=true;
		IO.log(GVars.logFile, 1, "Finished preprocessing of the reads", true);
		IO.writeToCommandLineBlockClose("FINISHED PRE-PROCESSING");
		return true;
	}
	
	
	public static void iterative5pTrimming(){
		
		if(GVars.species != null){
			
			Bowtie.resetMappingParameters();
			GVars.tempBowtieReportType = "-k";
			GVars.tempBowtieReportCount = "2";
			GVars.tempAlingType = "n";
			GVars.tempMM = 0;
			String tmpReads = GVars.output+File.separator+"reads.tmp";
			IO.copy(GVars.input, tmpReads, false);
			/**
			 * this will be the new input file
			 */
			String tempOut = GVars.output+File.separator+"reads_newIn.fa";
			
			for(int i = 0; i < GVars.iterative5pTrimming; i++){
				 
				AlignData dataHR = Bowtie.genomeAlignBowtie(tmpReads, GVars.output+File.separator+"genome_tmp.txt",true, true, false, true);
				if(dataHR != null){
					int r = i+1;
					IO.writeToCommandLineL1("Mapped in "+r+" round ("+i+" nucs removed from 5') of iterative 5' trimming: "+dataHR.rc);
					IO.log(GVars.logFile, 2, "Mapped in "+r+" round ("+i+" nucs removed from 5') of iterative 5' trimming: "+dataHR.rc, true);
					Bowtie.makeReadsFasta(dataHR.parsedFile, tempOut);
					removeIndexFasta(tmpReads, GVars.tmp, 1);
					new File(GVars.tmp).renameTo(new File(tmpReads));
				}
			}
			
			Bowtie.resetMappingParameters();
			// align using the user given parameters
			IO.copy(tmpReads, tempOut, true);

			boolean trimmed = GVars.adapterTrimmed;
			GVars.adapterTrimmed = true;
			boolean check =  readInput (tempOut, GVars.output+File.separator+"reads.fa",
					GVars.output+File.separator+"reads_withAdapter_dummy.fa",
					GVars.output+File.separator+"short_reads_dummy.txt","fasta");


			GVars.adapterTrimmed = trimmed;
			
			if(new File(GVars.output + File.separator+"reads.fa").exists() && reads > 0){
				GVars.input = GVars.output + File.separator+"reads.fa";
				IO.copy(GVars.output+File.separator+"reads.fa", GVars.output+File.separator+"reads_orig.fa", false);
				Util.sortSRNAbenchFormat(GVars.output+File.separator+"reads_orig.fa");
				Util.sortSRNAbenchFormat(GVars.output+File.separator+"reads.fa");				
			}
			else{
				IO.log(GVars.logFile, 4, "No reads found for input file: "+new File(GVars.input).getName()+". Will quit now.", true);
				IO.warning("No reads found for input file: "+new File(GVars.input).getName()+". Will quit now.");
				System.exit(0);
			}
			
			printPreProcResults();
			
		}
		else{
			IO.warning("Iterative 5' trimming can be only done when an assembly is given (species=)");
			IO.log(GVars.logFile, 4, "Iterative 5' trimming can be only done when an assembly is given (species=)", true);
		}
	}
	

	/**
	 * reads the input file and performs the preprocessing of the reads
	 * @param config
	 * @return
	 */
	public static boolean input( ){


		if(GVars.removeBarcode > 0){

			if(GVars.input.endsWith("fq") ||GVars.input.endsWith("fastq") || GVars.input.endsWith("FASTQ") || GVars.input.endsWith("fastQ") || 
					GVars.input.endsWith("fq.gz") ||GVars.input.endsWith("fastq.gz") 
					|| GVars.input.endsWith("FASTQ.gz") || GVars.input.endsWith("fastQ.gz")){

				removeIndex(GVars.input,GVars.input+".bcRem", GVars.removeBarcode );
				GVars.input = GVars.input+".bcRem";
				GVars.origInput = GVars.input+".bcRem";
			}
			else{
				IO.warning("The barcode trimming can be only applied to fastq files!");
				IO.log(GVars.logFile, 4, "The barcode trimming can be only applied to fastq files!", true);
				System.exit(0);
			}
		}
		

		/**
		 * will generate a read.fa file --> collaps all "redundant" reads
		 */
		if((GVars.inputType != null && GVars.inputType.equals("fastq")) || ( GVars.input.endsWith("fastq") 
				|| GVars.input.endsWith("FASTQ") || GVars.input.endsWith("fastQ")  || GVars.input.endsWith("fastq.gz") 
				|| GVars.input.endsWith("fq.gz")  || GVars.input.endsWith("fq") || GVars.input.endsWith(".bcRem"))){

			boolean check =  readInput (GVars.input, GVars.output+File.separator+"reads.fa",
					GVars.output+File.separator+"reads_withAdapter.fa",
					GVars.output+File.separator+"short_reads.txt","fastq");
			GVars.inputType = "fastq";
			
		}

		else if((GVars.inputType != null && GVars.inputType.equals("fasta")) || (GVars.input.endsWith("fa")
				|| GVars.input.endsWith("fasta") || GVars.input.endsWith("fas") || GVars.input.endsWith("fa.gz") || GVars.input.endsWith("fasta.gz")  || 
				GVars.input.endsWith("csfasta")  || GVars.input.endsWith("csfasta.gz")) ){
			
			boolean check =  readInput (GVars.input, GVars.output+File.separator+"reads.fa",
					GVars.output+File.separator+"reads_withAdapter.fa",
					GVars.output+File.separator+"short_reads.txt","fasta");
			GVars.inputType = "fasta";

		}
		
		else if((GVars.inputType != null && GVars.inputType.equals("sRNAbench")) || (GVars.input.endsWith("miR") || GVars.input.endsWith("detected") || GVars.input.endsWith("reads.fa"))){

			String[] in = GVars.input.split(":");
			new File(GVars.output+File.separator+"reads.fa").delete();
			for(String f : in){
				IO.copy(f, GVars.output+File.separator+"reads.fa",true);
			}
			GVars.inputType = "sRNAbench";
		}
		
		else if((GVars.inputType != null && GVars.inputType.equals("bowtie"))  || GVars.input.endsWith("parsed") || GVars.input.endsWith("bowtieOut")){
			
			GVars.adapterTrimmed = true;
			boolean check =  readInput (GVars.input, GVars.output+File.separator+"reads.fa",
					GVars.output+File.separator+"reads_withAdapter.fa",
					GVars.output+File.separator+"short_reads.txt","bowtieOut");
			GVars.inputType = "bowtie";

		}
		
		else{
			boolean check =  readInput (GVars.input, GVars.output+File.separator+"reads.fa",
					GVars.output+File.separator+"reads_withAdapter.fa",
					GVars.output+File.separator+"short_reads.txt","rc");
			GVars.inputType = "rc";

		}

		//////////////////////////////////////////////////////////////
		////// copy reads.fa to original 
		
		if(new File(GVars.output + File.separator+"reads.fa").exists() && reads > 0){
			GVars.input = GVars.output + File.separator+"reads.fa";
			IO.copy(GVars.output+File.separator+"reads.fa", GVars.output+File.separator+"reads_orig.fa", false);
			Util.sortSRNAbenchFormat(GVars.output+File.separator+"reads_orig.fa");
		}
		else{
			IO.log(GVars.logFile, 4, "No reads found for input file: "+new File(GVars.input).getName()+". Will quit now.", true);
			IO.warning("No reads found for input file: "+new File(GVars.input).getName()+". Will quit now.");
			System.exit(0);
		}
		
		/////////////////////////////////////////////////////////
		/// apply removal of 5' bases from adapter cleaned input
		if(GVars.remove3pBases > 0){
			reads = 0;
			readsUnique = 0;
			IO.log(GVars.logFile, 1, "Will remove the last "+GVars.remove3pBases+" nucleotides from 3' end of adapter trimmed reads", true);
			IO.writeToCommandLineL2("Will remove the last "+GVars.remove3pBases+" nucleotides from 3' end of adapter trimmed reads");
			remove3pBases(GVars.input, GVars.output+File.separator+"short_reads.txt",GVars.remove3pBases);
			IO.copy(GVars.output+File.separator+"reads.fa", GVars.output+File.separator+"reads_orig.fa", false);
			Util.sortSRNAbenchFormat(GVars.output+File.separator+"reads_orig.fa");
		}
		
		/// print out the results of the preprocessing
		printPreProcResults();
		return true;
	}
	
	
	/**
	 * This function is applied to a reads.fa file (sRNAbench internal format)<br>
	 * It removes the last bases from the 3' end of the adapter trimmed reads 
	 * @param file
	 * @param shortReadsFile
	 * @param bases
	 */
	public static void remove3pBases(String file, String shortReadsFile, int bases){
		
		Map<String,int[]> count = new Hashtable<String,int[]>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			BufferedWriter writer = new BufferedWriter (new FileWriter(GVars.tmp));
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("#");
				int rc = Integer.parseInt(f[1].trim());
				String seq = reader.readLine();
				if(seq.length() - bases > 0){
					String subSeq = seq.substring(0, seq.length()-bases);
					Util.addIntMap(count, subSeq, rc);
				}
			}

			reader.close();
			writer.close();
			writeOutReadsFile(count, file, shortReadsFile, true);
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Input file not found (Preproc.remove3pBases)"+new File(GVars.input).getName(), true);
			IO.warning("Input file not found (Preproc.remove3pBases)"+new File(GVars.input).getName());
			System.exit(1);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "tmp file "+new File(GVars.tmp).getName()+" (Preproc.remove3pBases) cannot be written", true);
			IO.warning("tmp file "+new File(GVars.tmp).getName()+" (Preproc.remove3pBases) cannot be written");
			System.exit(1);
		}

		
	}
	
	
	/**
	 * Allows to remove the first bases of an fastq file (barcode)
	 * @param infile
	 * @param outfile
	 * @param firstBases
	 */
	public static void removeIndex(String infile, String outfile, int firstBases){

	

			String[] files = infile.split(":");
			BufferedReader reader;


			try {

				BufferedWriter writer = new BufferedWriter (new FileWriter(outfile));
				for(String file1 : files){

					System.out.println();
					IO.writeToCommandLineL2("Reading file: "+file1);
					if(file1.endsWith("gz")){
						reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
					}
					else{
						reader = new BufferedReader (new FileReader(file1));
					}
					String line = null;

					while((line = reader.readLine()) != null){



						String seq = reader.readLine();
						String dummy = reader.readLine();
						String qual = reader.readLine();


						writer.write(line+"\n");
						writer.write(seq.substring(firstBases, seq.length())+"\n");
						writer.write(dummy+"\n");
						writer.write(qual.substring(firstBases,  qual.length())+"\n");

					}				

					reader.close();
				}
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
	
	
	/**
	 * Allows to remove the first bases of an fastq file (barcode)
	 * @param infile
	 * @param outfile
	 * @param firstBases
	 */
	public static void removeIndexFasta(String infile, String outfile, int firstBases){

	

			String[] files = infile.split(":");
			BufferedReader reader;


			try {

				BufferedWriter writer = new BufferedWriter (new FileWriter(outfile));
				for(String file1 : files){

					System.out.println();
					IO.writeToCommandLineL2("Reading file: "+file1);
					if(file1.endsWith("gz")){
						reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
					}
					else{
						reader = new BufferedReader (new FileReader(file1));
					}
					String line = null;

					while((line = reader.readLine()) != null){


						if(line.startsWith(">")){
							String seq = reader.readLine();
							writer.write(line+"\n");
							writer.write(seq.substring(firstBases, seq.length())+"\n");
						}

					}				

					reader.close();
				}
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
	
	
	/**
	 * gets the position of the first base that has a quality below minQ
	 * @param qual
	 * @param minQ
	 * @param rest
	 * @return
	 */
	public static int getNumberOfnucsBelowMinQ(String qual, int minQ, int phredRest){

		char[] q = qual.toCharArray();
		int count=0;
		for(int i = 0; i < q.length; i++){
			int val = q[i];
			if(val-phredRest < minQ)
				count++;

		}
		return count;
	}
	/**
	 * gets the position of the first base that has a quality below minQ
	 * @param qual
	 * @param minQ
	 * @param rest
	 * @return
	 */
	public static double getMeanQ(String qual, int phredRest){

		char[] q = qual.toCharArray();
		int count=0;
		for(int i = 0; i < q.length; i++){
			int val = q[i];
			count+= (val - phredRest);
		}
		return (double)count/(double)q.length;
	}
	
	
	
	/**
	 * gives back true if the quality criterion is fulfilled
	 * @param qual
	 * @param rc
	 * @return
	 */
	public static boolean checkQuality(String qual, int rc){
		
		
		if(qual != null && GVars.qualityType.equals("min")){
			int  nr =  getNumberOfnucsBelowMinQ(qual, GVars.minQ, GVars.phred);
			if(nr > GVars.maxQfailure){
				readsQRCfiltered+= rc;
				return false;
			}
		}
		else if(qual != null && GVars.qualityType.equals("mean")){
			double mean = getMeanQ(qual, GVars.phred);
			if(mean < GVars.minQ){
				readsQRCfiltered+= rc;
				return false;
			}
		}
		return true;
	}
	
	/**
	 * takes a read and tries to detect the adapter
	 * @param read
	 * @param map
	 * @return: true --> adapter was found; false --> adapter was not found
	 * 
	 * @uses: GVars.recursiveAdapterTrimming; GVars.holdNonAdapter
	 */
	public static boolean trimAdapter(String read, String qual, Map<String,int[]> map, int rc){

		// input is adapter trimmed --> add read 
		if(GVars.adapterTrimmed){

			if(checkQuality(qual, rc)){
				Util.addIntMap(map, new String(read), rc);
				return true;
			}
			else
				return false;
		}

		int pos = Align.alignAdapter(read, GVars.adapter, GVars.adapterStart, GVars.adapterMM, GVars.adapterMinLength);
		// if adapter is detected --> add to map and return true
		if(pos >= 0){
			
			String seq = read.substring(0, pos);
			String nQual = null;
			if(qual != null)
				nQual = qual.substring(0, pos);
			
			if(checkQuality(nQual,rc)){
				Util.addIntMap(map, seq, rc);
				readsAdapterFound++;
				return true;
			}
			else
				return false;
			
		}
		
		
		if(GVars.recursiveAdapterTrimming){

			pos = recursiveAdapterTrimming(read);
			//// in recursive adapter trimming, all reads are keept in the analysis

			if(pos >= 0){

				String seq = read.substring(0, pos);
				
				String nQual = null;
				if(qual != null)
					nQual = qual.substring(0, pos);
				
				if(checkQuality(nQual,rc)){
					Util.addIntMap(map, seq, rc);
					readsAdapterFound++;
					return true;
				}
				else
					return false;
				
			}
			else{
				
				
				if(checkQuality(qual,rc)){
					Util.addIntMap(map, new String (read), rc);
					return true;
				}
				else
					return false;
				
			}
		}
		
		readsAdapterNotFound++;

		// if the sequences where the adapter has not been found are used 
		if(GVars.holdNonAdapter){
			
			if(checkQuality(qual,rc)){
				Util.addIntMap(map, new String (read), rc);
				return true;
			}
			else
				return false;
			
		}
		return false;
	}
	
	
	
	/**
	 * performs recursive adapter trimming
	 * @param read
	 * @return: -1 --> adapter was not found; >= 0 --> adapter was found at that position (in the read string)
	 */
	public static int recursiveAdapterTrimming(String read){
		
		int pos = -1;
		for(int i = read.length() - GVars.adapterMinLength +1 ; i < read.length();i++){
			pos = Align.alignAdapter(read, GVars.adapter, i, 0, read.length() - i );

			if(pos >= 0)
				return pos;
		}
		return -1;
	}
	
	
	/**
	 * the input format is fasta. Will read one more line. 
	 * @param line
	 * @param reader
	 * @param writerNoAdapter
	 * @param file
	 * @return
	 */
	public static String[] isFasta(String line, BufferedReader reader, BufferedWriter writerNoAdapter, String file){
		
		String seq = null;
		int rc = 0;
		try {
			seq = reader.readLine();
			if(seq == null){
				IO.warning("Read Sequence expected after line: "+line+" but not found in file "+new File(file).getName());
				IO.log(GVars.logFile, 3, "Read Sequence expected after line: "+line+" but not found in file "+new File(file).getName(), true);
				return null;
			}
			seq = seq.trim();

			String[] f = line.split("\\s+");
			String[] cc = f[0].split(GVars.sep);	
			
			if(f.length >= 2){
				rc = Integer.parseInt(f[1]);
			}
			else{
				if(cc.length >= 2)
					rc = Integer.parseInt(cc[1]);
				else
					rc = 1;
			}
			
			if(f.length < 2 && cc.length < 2){
//				IO.warning("Error in input file. Expected was a fasta file with "+GVars.sep+" (or space) as separator like ID"+GVars.sep+"RC being RC the read count (a number).");
//				IO.log(GVars.logFile, 4, "Error in input file. Expected was a fasta file with "+GVars.sep+" (or space) as separator like ID"+GVars.sep+"RC being RC the read count (a number).", true);
//				return null;
				IO.log(GVars.logFile, 4, "Warning for fasta input file. Expected was a fasta file with "+GVars.sep+" (or space) as separator like ID"+GVars.sep+"RC being RC the read count (a number).", true);
			}

			if(GVars.writeNonAdapter && writerNoAdapter != null){
				writerNoAdapter.write(line+"\n");
				writerNoAdapter.write(seq+"\n");
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Error in file (unexpected end of file?): ", true);
			IO.warning("Error in input file (unexpected end of file?): ");
			return null;
//			e.printStackTrace();
		}
		catch (Exception e){
			IO.log(GVars.logFile, 4, "General exception error in line: "+line, true);
			IO.warning( "General exception error in line: "+line);
			return null;
			
		}


		String[] back = new String[3];
		back[0] = seq;
		back[1] = rc+"";
		return back;

	}
	
	/**
	 * The input is in read/count format
	 * @param line
	 * @param reader
	 * @param writerNoAdapter
	 * @param file
	 * @return
	 */
	public static String[] isReadCount(String line, BufferedReader reader, BufferedWriter writerNoAdapter, String file){
		
		int rc = 0;
		String seq = null;
		try {

			String[] f = line.split("\t");
			if(f.length >= 2){
				rc = Integer.parseInt(f[1]);
				seq = f[0];
			}


			if(GVars.writeNonAdapter && writerNoAdapter != null){
				writerNoAdapter.write(line+"\n");
			}

		}catch(NumberFormatException e){
			IO.log(GVars.logFile, 3, "No integer found in second column of read/count format file (probably it is the header).", true);
			IO.warning("No integer found in second column of read/count format file (probably it is the header).");
			return "NA:0".split(":");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Error in file (unexpected end of file?): ", true);
			IO.warning("Error in input file (unexpected end of file?): ");
			return null;
//			e.printStackTrace();
		}
		catch (Exception e){
			IO.log(GVars.logFile, 4, "General exception error in line: "+line, true);
			IO.warning( "General exception error in line: "+line);
			return null;
			
		}


		String[] back = new String[3];
		back[0] = seq;
		back[1] = rc+"";
		return back;

	}
	
	/**
	 * the input file is a bowtie output file
	 * @param line
	 * @param file
	 * @return: has an additional back value [2] --> the ID of the read
	 */
	public static String[] isBowtieOut(String line, String file){
		
		int rc = 0;
		String seq = null;
		String[] f = null;
		try {

			f = line.split("\t");
			
			seq = f[4];
			if(f[1].equals("-"))
				seq = SeqUtil.getReverseComplementarySequence(seq);
			
			String[] f1 = f[0].split("#");
			rc = 1;
			if(f1.length >= 2)
				rc = Integer.parseInt(f1[1]);



		}catch(NumberFormatException e){
			return null;
		} 

		String[] back = new String[3];
		back[0] = seq;
		back[1] = rc+"";
		back[2] = f[0];
		return back;

	}

	
	/**
	 * the input is fastq. will read 3 more lines 
	 * @param line: the first line of the 4 line fastq block
	 * @param reader: the BufferedReader object
	 * @param writerNoAdapter --> write out non-adapter trimmed reads
	 * @param file: the input file --> only for loggin in case of an error
	 * @return; a String[]: [0] --> the read sequence, [0] --> the read count (always 1 in case of fastq)
	 */
	public static String[] isFastq(String line, BufferedReader reader, BufferedWriter writerNoAdapter, String file){
		
		String seq = null;
		String qual = null;
		try {
			seq = reader.readLine();
			if(seq == null){
				IO.warning("Read Sequence expected after line: "+line+" but not found in file "+new File(file).getName());
				IO.log(GVars.logFile, 4, "Read Sequence expected after line: "+line+" but not found in file "+new File(file).getName(), true);
				return null;
			}
			seq = seq.trim();

			String idDumm = reader.readLine();
			qual = reader.readLine();
			if(qual == null){
				IO.warning("Input file "+new File(file).getName()+" has not a multiple of 4!");
				IO.log(GVars.logFile, 4, "Input file "+new File(file).getName()+" has not a multiple of 4!", true);
				return null;
			}
			
			
			if(GVars.writeNonAdapter && writerNoAdapter != null){
				writerNoAdapter.write(line+"\n");
				writerNoAdapter.write(seq+"\n");
				writerNoAdapter.write(idDumm+"\n");
				writerNoAdapter.write(qual+"\n");
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Error in file (unexpected end of file?): ", true);
			IO.warning("Error in input file (unexpected end of file?): ");
			return null;
//			e.printStackTrace();
		}
		catch (Exception e){
			IO.log(GVars.logFile, 4, "General exception error in line: "+line, true);
			IO.warning( "General exception error in line: "+line);
			return null;
			
		}

		String[] back = new String[3];
		back[0] = seq;
		back[1] = 1+"";
		back[2] = qual;
		return back;
	}
	
	
	/**
	 * Trim the adapter from a fastq file
	 * @param fastq
	 * @param output
	 * @return
	 */
	public static int trimmFastq(String fastq, String output, double ratio){

		String[] files = fastq.split(":");
		int c = 0;
		int f = 0;
		BufferedReader reader;

		try {

			BufferedWriter writer = new BufferedWriter(new FileWriter(output));

			for(String file1 : files){

				System.out.println();
				IO.writeToCommandLineL2("Reading file: "+file1);

				if(file1.endsWith("gz")){
					reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
				}
				else{
					reader = new BufferedReader (new FileReader(file1));
				}
				try{

					String line = null;

					while((line = reader.readLine()) != null){

						String seq = reader.readLine();
						String dumm = reader.readLine();
						String qual = reader.readLine();
						if(seq != null && qual != null){
							double ratioS = SeqUtil.ratioOfMostFrequentBase(seq);
							if(ratioS >= ratio)
								continue;

							int pos = Align.alignAdapter(seq, GVars.adapter, GVars.adapterStart, GVars.adapterMM, GVars.adapterMinLength);
							if(pos >= 0){
								seq = seq.substring(0, pos);
								qual = qual.substring(0, pos);
								f++;
							}
							if(seq.length() >= GVars.minReadLength){
								writer.write("@"+c+"\n");
								writer.write(seq+"\n");
								writer.write("+"+c+"\n");
								writer.write(qual+"\n");
								c++;
							}
		
						}
					}
					writer.close();
				}

				catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		} catch(IOException e){

		}
		IO.writeToCommandLineL1("Found "+f+" reads with adapter. "+c+" reads in output file. ");
		return c;
	}

	/**
	 * 
	 * @param file
	 * @param readOutFile
	 * @param noAdapterFile
	 * @param lengthFilteredOut
	 * @param type: the input file type. can be 'fastq', 'fasta', 'rc'
	 * @return
	 * 
	 * @uses: all adapter parameters from GVars
	 */
	public static boolean readInput (String file, String readsFile, String noAdapterFile, String shortReadsFile, String type){


//		Map<Thread, TrimmAdapter> threads = new Hashtable<Thread,TrimmAdapter>();
		String[] files = file.split(":");

		Set<String> detected = new HashSet<String>();
		Map<String,int[]> map = new Hashtable<String,int[]>();
		BufferedReader reader;

		try {

			BufferedWriter writerNoAdapter = new BufferedWriter(new FileWriter(noAdapterFile));

			for(String file1 : files){

				System.out.println();
				IO.writeToCommandLineL2("Reading file: "+file1);
				
				if(file1.endsWith("gz")){
					reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
				}
				else{
					reader = new BufferedReader (new FileReader(file1));
				}
				try{
					
					String line = null;

					while((line = reader.readLine()) != null){

						String[] data = null;
						if(type.equals("fastq")){
							data = isFastq(line, reader, writerNoAdapter,  file1);
						}
						else if(type.equals("fasta")){
							data = isFasta(line, reader, writerNoAdapter, file1);
						}
						else if(type.equals("rc")){
							data = isReadCount(line, reader, writerNoAdapter, file1);
						}
						else if(type.equals("bowtieOut")){
							data = isBowtieOut(line, file1);
							if(data == null)
								break;
							if(detected.contains(data[2])){
								continue;
							}
							else{
								detected.add(data[2]);
							}
						}
						
						
						if(data == null){
							break;
						}
						else{
							if(Integer.parseInt(data[1]) == 0){
								break;
							}
						}
						// number of input reads
						readsRaw += Integer.parseInt(data[1]);
						if(readsMaxLengthInput < data[0].length())
							readsMaxLengthInput = data[0].length();

//						System.out.println(map.size());
//						if(threads.size() < GVars.p){
//							
//							TrimmAdapter ta = new TrimmAdapter(new String(data[0]),Integer.parseInt(data[1]));
//							Thread t = new Thread(ta);					
//							t.start();
//							threads.put(t, ta);	
//						}
//						else{
//							
//							while(threads.size() >= GVars.p){
//								checkTreads(threads,map);
//							}
//							// add the current read once the number of 
//							TrimmAdapter ta = new TrimmAdapter(new String(data[0]),Integer.parseInt(data[1]));
//							Thread t = new Thread(ta);
//							t.start();
//							threads.put(t, ta);	
//						}

						
						if(GVars.qualityType == null){
							trimAdapter(data[0], null, map, Integer.parseInt(data[1]));							
						}
						else{
							trimAdapter(data[0], data[2], map, Integer.parseInt(data[1]));							
						
						}
///////////////////////////////////////////////////////////////

					}

					reader.close();
				} 
				catch (FileNotFoundException e){
					System.out.println(file+" not found");
					System.exit(1);
				}
			}
			writerNoAdapter.close();
			if(!(GVars.writeNonAdapter)){
				new File(noAdapterFile).delete();
			}



		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			IO.warning("File not found or error in reading file: "+new File(file).getName()+". Will quit now. ");
			IO.log(GVars.logFile, 3, "File not found or error in reading file: "+new File(file).getName()+". Will quit now. ", true);
			System.exit(1);
		}

//		while(threads.size() > 0){
//			checkTreads(threads,map);
//		}
		writeOutReadsFile(map, readsFile, shortReadsFile, false);
		return true;
	}

	
//	/**
//	 * This functions: i) checkes if the threads have finished already, ii) adds the result to the map, iii) removes the finished Threads from the Map
//	 * @param threadMap
//	 * @param map
//	 */
//	private static void checkTreads(Map<Thread,TrimmAdapter> threadMap, Map<String,int[]> map ){
//		
//		for(Iterator<Map.Entry<Thread,TrimmAdapter>> it = threadMap.entrySet().iterator(); it.hasNext();){
//			Map.Entry<Thread,TrimmAdapter> entry = it.next();
//			if(entry.getKey().isAlive()){
//				
//			}
//			else{
//				if(entry.getValue().add){
//					Util.addIntMap(map, entry.getValue().adapterTrimmedRead , entry.getValue().rc);
//					readsAdapterFound += entry.getValue().readsAdapterFound; 
//					readsAdapterNotFound += entry.getValue().readsAdapterNotFound; 
//				}
//				it.remove();
//			}
//		}
//		
//	}
	
	/**
	 * Writes out the internal sRNAbench fasta format with collapsed reads
	 * @param map: keys:read sequence; values: int[0]=read count 
	 * @param readsFile: the reads.fa file
	 * @param shortReadsFile: the file with the 'to short reads'
	 * @param append2ShortReads: true --> append to the shortsReadsFile (if exists) 
	 */
	private static void writeOutReadsFile(Map<String,int[]> map, String readsFile, String shortReadsFile, boolean append2ShortReads){
		
		int c= 1;
		try {
			BufferedWriter readWriter  = new BufferedWriter(new FileWriter(readsFile));
			BufferedWriter filteredWriter = new BufferedWriter(new FileWriter(shortReadsFile,append2ShortReads));
			filteredWriter.write("read\tcount\n");

			for(String seq : map.keySet()){	


				if(!(seq.contains("N") || seq.contains(".")) && 
						map.get(seq)[0] >= GVars.minRC && seq.length() >= GVars.minReadLength && seq.length() <= GVars.maxReadLength){
					
					readWriter.write(">"+c+"#"+map.get(seq)[0]+"\n");
					readWriter.write(seq+"\n");
					c++;
					
					reads += map.get(seq)[0];
					readsUnique++;
					
					// set the maximum length
					if(readsMaxLengthAnalysis < seq.length())
						readsMaxLengthAnalysis = seq.length();
				}
				else if(!(seq.contains("N") || seq.contains(".")) && 
						map.get(seq)[0] >= GVars.minRC && seq.length() < GVars.minReadLength){
					readsLengthFilteredMin += map.get(seq)[0];
					filteredWriter.write(seq+"\t"+map.get(seq)[0]+"\n");
				}
				else if(!(seq.contains("N") || seq.contains(".")) && 
						map.get(seq)[0] >= GVars.minRC && seq.length() > GVars.maxReadLength){
					readsLengthFilteredMax +=map.get(seq)[0];
				}
				else{
					readsQRCfiltered += map.get(seq)[0];
				}
			}
			//			back[1]=back[0]-(back[2]+back[3]);
			readWriter.close();
			filteredWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	////// SOLiD specific functions 

	/**
	 * 1) preprocessing of the input reads (adapter trimming and collapsing) <br>
	 * 2) the reads in colour space are aligned to the genome <br>
	 * 3) extract the reads in nucleotide space and re-collaps them (unique reads in colour space need not to be unique in nuc. space). <br>
	 */
	private static void prepareSolid (){
		input();
		deleteInfo();
		if( GVars.species != null){

			if(GVars.doAlignment){
				IO.writeToCommandLineBlockOpen("PREPARING SOLID");
				IO.log(GVars.logFile, 1, "PREPARING SOLID", true);
				AlignData ad = Bowtie.genomeAlignBowtie(GVars.input,GVars.output+File.separator+"genome.txt", true, true, false, true);
				// if solid --> extract sequences in nuc. space & realign
				if(ad != null){
					fixParsedFileSolid(GVars.output+File.separator+"genome.parsed", GVars.output+File.separator+"readsC.fa");
					Util.sortSRNAbenchFormat(GVars.output+File.separator+"readsC.fa");
					GVars.colorFlag="";
					GVars.colorIndex = "";
					GVars.solid = false;
					GVars.input = GVars.output+File.separator+"readsC.fa";
//					GVars.holdNonAdapter =true;
//					GVars.adapter = null;
					GVars.inputType = "fasta";
					GVars.adapterTrimmed = true;
					//					check = Align.genome();		

				}
				else{
					IO.log(GVars.logFile, 4, "No reads have been mapped to the genome.", true);
					IO.warning("No reads have been mapped to the genome (Preproc.prepareSolid()). Will quit now.");
					System.exit(1);
				}
				IO.writeToCommandLineBlockClose("FINISHED PREPARING SOLID");
			}
		}
	}
	
	/**
	 * unique reads in colour space need not to be unique in nuc. space! This function regroups them!
	 * @param parsed: the bowtie output file
	 * @param: the 'new' input file for sRNAbench (in sRNAbench format!)
	 */
	public static void fixParsedFileSolid(String parsed, String out){
		
		Set<String> detected = new HashSet<String>();
		Map<String,int[]> reads = new Hashtable<String,int[]>();
		try {
			BufferedReader reader = new BufferedReader ( new FileReader(parsed));
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				int rc = Integer.parseInt(f[0].split("#")[1]);
				String read = f[4];
				if(f[1].equals("-"))
					read = SeqUtil.getReverseComplementarySequence(read);
				
				if(detected.contains(f[0])){
					
				}
				else{
					detected.add(f[0]);
					if(reads.containsKey(read)){
						reads.get(read)[0] += rc;
					}
					else{
						int[] t = new int[1];
						t[0] = rc;
						reads.put(read, t);
					}
					
				}
				
				
			}
			reader.close();
			BufferedWriter writer = new BufferedWriter( new FileWriter(out));
			int c = 1;
			for(String key : reads.keySet()){
				writer.write(">"+c+"#"+reads.get(key)[0]+"\n");
				writer.write(key+"\n");
				c++;
			}
			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

	///////////  	END PREPARE SOLID
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * 
	 * @param file
	 * @param minAdapterLength
	 * @param startPos
	 * @return: the adapter sequence
	 */
	public static String  getMostLikelyAdapter (String file,  int minAdapterLength, int startPos, String type){

		String[] files = file.split(":");

		Map<String,int[]> map = new Hashtable<String,int[]>();
		BufferedReader reader;

		try {
			for(String file1 : files){
				if(file.endsWith("gz")){
					reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
				}
				else{
					reader = new BufferedReader (new FileReader(file1));
				}
				String line = null;

				while((line = reader.readLine()) != null){

					String[] data = null;
					if(type.equals("fastq")){
						data = isFastq(line, reader, null,  file1);
					}
					else if(type.equals("fasta")){
						data = isFasta(line, reader, null, file1);
					}
					else if(type.equals("rc")){
						data = isReadCount(line, reader, null, file1);
					}
					
					
					if(data == null){
						break;
					}


					for(int i = startPos; i <= data[0].length() - minAdapterLength; i++){
						String ss = data[0].substring(i, i+minAdapterLength);
						Util.addIntMap(map, ss, 1);

					}
				}
				reader.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		String back = "A";
		int highest = 0;
		for(String key : map.keySet()){
			if(SeqUtil.ratioOfMostFrequentBase(key) < 0.7 && !(key.contains("N"))){
				if(map.get(key)[0] > highest){
					back = key;
					highest = map.get(key)[0];
				}
			}
		}
		return new String(back);
	}



	
	/**
	 * Prints the preprocessing results to screen
	 */
	private static void printPreProcResults(){
		
		DecimalFormat df = new DecimalFormat("#.##");

		IO.writeToCommandLineL1("Result of pre-processing");
		IO.writeToCommandLineL2("No. raw input reads: "+readsRaw); 
		IO.log(GVars.logFile, 2, "No. raw input reads: "+readsRaw, true);

		double ratReads = 100d*(double)reads / (double)readsRaw;
		double ratReadsU = 100d*(double)readsUnique / (double)readsRaw;

		IO.writeToCommandLineL2("Number of reads in analysis: "+reads+" ("+df.format(ratReads)+"% of raw input reads)");
		IO.log(GVars.logFile, 2, "Number of reads in analysis: "+reads+" ("+df.format(ratReads)+"% of raw input reads)", true);

		
		IO.writeToCommandLineL2("Number of unique reads in analysis: "+readsUnique+" ("+df.format(ratReadsU)+"% of raw input reads)");		
		IO.log(GVars.logFile, 2, "Number of unique reads in analysis: "+readsUnique+" ("+df.format(ratReadsU)+"% of raw input reads)", true);


		if(GVars.adapterTrimmed){
			IO.writeToCommandLineL2("Input was adapter trimmed");

		}
		else{

			IO.writeToCommandLineL2("No. input reads where the adapter was not found: "+readsAdapterNotFound);
			IO.log(GVars.logFile, 2, "No. input reads where the adapter was not found: "+readsAdapterNotFound, true);
		
			double adapterPerc = 100d*(double)readsAdapterFound/ (double)readsRaw;

			IO.writeToCommandLineL2("No. adapter trimmed input reads: "+readsAdapterFound+" ("+df.format(adapterPerc)+" % of raw input reads)");
			IO.log(GVars.logFile, 2, "No. adapter trimmed input reads: "+readsAdapterFound+" ("+df.format(adapterPerc)+" % of raw input reads)", true);
		}

		IO.writeToCommandLineL2("No. length filtered input reads (min. Length): "+readsLengthFilteredMin);	 
		IO.log(GVars.logFile, 2, "No. length filtered input reads (min. Length): "+readsLengthFilteredMin, true);

		IO.writeToCommandLineL2("No. length filtered input reads (max. Length): "+readsLengthFilteredMax);	 
		IO.log(GVars.logFile, 2, "No. length filtered input reads (max. Length): "+readsLengthFilteredMax, true);

		IO.writeToCommandLineL2("Max. read length in input file: "+readsMaxLengthInput);	
		IO.log(GVars.logFile, 2, "Max. read length in input file: "+readsMaxLengthInput, true);

		IO.writeToCommandLineL2("Max. read length in analysis: "+readsMaxLengthAnalysis);	
		IO.log(GVars.logFile, 2, "Max. read length in analysis: "+readsMaxLengthAnalysis, true);

		IO.writeToCommandLineL2("Filtered reads (low quality or low read count): "+readsQRCfiltered);	
		IO.log(GVars.logFile, 2, "Filtered reads (low quality or low read count): "+readsQRCfiltered, true);
		
	}


	/**
	 * Aligns the input reads to fasta/bowtie-index libraries and removes the mapped reads from the analysis <br>
	 * This function will use -k 1 for the bowtie mapping as the read will be removed when it maps to one reference sequence. 
	 * 
	 * @param filterLibs: a list of libraries to be filtered
	 * @param readsFile: the file with the input reads
	 * @param overwrite: overwrite the bowtie index: true = overwrite; false = do not overwrite
	 */
	public static void applyFilter(List<String> filterLibs, String readsFile, String outPath, boolean overwrite){

		if(filterLibs != null){
			
			IO.log(GVars.logFile, 1, "Apply filter libraries", true);
			/////////////// only testing - remove later!!!!!!
			int[] c = Stat.getCountsFastaFile(readsFile, "#");
//			System.out.println(reads+" "+c[1]);
//			System.out.println(readsUnique+" "+c[0]);
			/////////////////


			for(String library : filterLibs){
				String[] index = Bowtie.getIndexAndName(library, overwrite);

				Bowtie.setMappingParamters(library);
				GVars.tempBowtieReportType = "-k";
				GVars.tempBowtieReportCount = "2";

				String parsedFile = Bowtie.mapToIndex(GVars.output,index[0], index[1], readsFile, GVars.colorFlag, GVars.colorIndex, true);

				if(parsedFile != null){
					AlignData ad = Stat.profileFromBowtieOut(parsedFile);
					ad.library = index[1];
					Set<String> set = Bowtie.getBowtieReferences(parsedFile, 0);
					Util.removeReadsFasta(readsFile, set);

				}


			}

			int[] c1 = Stat.getCountsFastaFile(readsFile, "#");
			IO.writeToCommandLineL1("After filtering "+c1[0]+" unique reads and "+c1[1]+" read count");	
		}

	}
	
	/**
	 * Overwrites the previous results 
	 */
	private static void deleteInfo(){

		readsRaw = 0;
		readsAdapterFound = 0;
		readsAdapterNotFound = 0;
		readsLengthFilteredMin = 0;
		readsLengthFilteredMax = 0;
		reads = 0;
		readsUnique = 0;
		readsQRCfiltered = 0;
		readsMaxLengthInput = 0;

	}
}
