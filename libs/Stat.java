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
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;



public class Stat {

	
	
	/**
	 * takes a fasta file in sRNAbench/miRanalyzer format and gives back the unique read number [0] and the read count [1]
	 * @param file: the fasta file
	 * @param sep: the separator between ID and read count, for example >1#299 --> read with ID=1 has read count 299 and is separated by #
	 * @return [0] --> UR; [1] --> RC
	 */
	public static int[] getCountsFastaFile(String file, String sep ){


		int[] c = new int[2];
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(file));
			String line = null;

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){

					String[] f = line.replace(">", "").trim().split(sep);
					if(f.length >= 2){
						c[1] += Integer.parseInt(f[1]);
						c[0]++;

					}
					else{
						c[0]++;
						c[1]++;
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return c;

	}
	
	/**
	 * 
	 * @param input
	 * @param readFile
	 */
	public static void makeReadLengthStat(String input, String readFile){
		
		/**
		 * will generate a read.fa file --> collaps all "redundant" reads
		 */
		String fullReadLength = GVars.stat+File.separator+"readLengthFull.txt";
		String analysisReadLength = GVars.stat+File.separator+"readLengthAnalysis.txt";
		
		if(GVars.inputType.equals("fastq")){

			writeReadLength (input,fullReadLength,"fastq",Preproc.readsMaxLengthInput);
			
		}

		else if(GVars.inputType.equals("fasta")){
			
			writeReadLength (input,fullReadLength,"fasta",Preproc.readsMaxLengthInput);	

		}
		
		else if(GVars.inputType.equals("rc")){
			writeReadLength (input,fullReadLength,"rc",Preproc.readsMaxLengthInput);
		}
		
		
		makeReadLengthSRNAbenchFormat(readFile, analysisReadLength);

		
		if(GVars.graphics){
			graphics.getGraph(GVars.stat+File.separator+"readLengthAnalysis.txt", GVars.graphs+File.separator+"readLengthAnalysis.png", 1);
			graphics.getGraph(GVars.stat+File.separator+"readLengthFull.txt", GVars.graphs+File.separator+"readLengthFull.png", 1);
		}
		
		
		
	}
	
	
	public static void makeReadLengthSRNAbenchFormat(String input, String output){
		
		boolean o = GVars.adapterTrimmed;
		GVars.adapterTrimmed = true;
		String oSep = GVars.sep;
		GVars.sep = "#";
		
		writeReadLength(input, output, "fasta", Preproc.readsMaxLengthAnalysis);
		// reestablish the user parameters
		GVars.sep = oSep;
		GVars.adapterTrimmed = o;
	}
	
	/**
	 * Reads the original input file and applies recursive adapter trimming to get the read length distribution
	 * @param file: the input file
	 * @param outFile: the read length output file
	 * @param type: the file type: fastq, fasta, rc
	 * @param maxReadLength
	 */
	public static void writeReadLength (String file, String outFile, String type, int maxReadLength){

		String[] files = file.split(":");
		BufferedReader reader;
		Set<String> detected = new HashSet<String>();
		try {

			int[] ur = new int[maxReadLength+1];
			int[] rc = new int[maxReadLength+1];
			Map<String,int[]> map = new Hashtable<String,int[]>();
			for(String file1 : files){

				
				if(file1.endsWith("gz")){
					reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
				}
				else{
					reader = new BufferedReader (new FileReader(file1));
				}
				try{
					
					String line = null;

					while((line = reader.readLine()) != null){

						//					catch if seq is null!!!
						String[] data = null;
						if(type.equals("fastq")){
							data = Preproc.isFastq(line, reader, null,  file1);
						}
						else if(type.equals("fasta")){
							data = Preproc.isFasta(line, reader, null, file1);
						}
						else if(type.equals("rc")){
							data = Preproc.isReadCount(line, reader, null, file1);
						}
						else if(type.equals("bowtieOut")){
							data = Preproc.isBowtieOut(line, file1);
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
						// number of input reads
						if(GVars.qualityType == null){
							Preproc.trimAdapter(data[0], null, map, Integer.parseInt(data[1]));							
						}
						else{
							Preproc.trimAdapter(data[0], data[2], map, Integer.parseInt(data[1]));							
						
						}

					}

					reader.close();

				}
				catch (FileNotFoundException e){
					System.out.println(file+" not found");
					System.exit(1);
				}
			}
			
			for(String key : map.keySet()){
				rc[key.length()] += map.get(key)[0];
				ur[key.length()] ++;
			}
			writeReadLength(outFile, ur, rc);

		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * writes out the read length distribution file
	 * @param outFile
	 * @param ur: the number of unique reads as a function of length (the index)
	 * @param rc: the number of reads as a function of length (the index)
	 */
	private static void writeReadLength(String outFile, int[] ur, int[] rc ){
		
		BufferedWriter writer;
		try {
			writer = new BufferedWriter(new FileWriter(outFile));
			writer.write("Read Length (nt)\tUR\tPercentage_UR\tRC\tPercentage_RC\tRPM\n");
			int totUR = getSum(ur);
			int totRC = getSum(rc);
			
			for(int i = 0; i < ur.length; i++){

				double urt = 100d*(double)ur[i]/(double)totUR;
				double rdt = 100d*(double)rc[i]/(double)totRC;
				double rpm = (double)rc[i]/((double)totRC/1000000d);
				writer.write(i+"\t"+ur[i]+"\t"+urt+"\t"+rc[i]+"\t"+rdt+"\t"+rpm+"\n");
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 3, "Error when writing "+new File(outFile).getName()+" (Stat.writeReadLength)", true);
			IO.warning("Error when writing "+new File(outFile).getName()+" (Stat.writeReadLength)");
		}

	}
	
	
	/**
	 * gives back the sum of an int array
	 * @param array
	 * @return
	 */
	public static int getSum(int[] array){
		int back = 0;
		for(int i : array){
			back+=i;
		}
		return back;
	}


	
	/**
	 * get mapping counts as a function of reference name from bowtie output file
	 * @param bowtieFile
	 * @return String <-->double[]; the array holds: [0]-->UR+; [1]-->RC+;[2]-->multiple mapping adjusted read count(+) <br>
	 *  [3]-->UR (-); [4]-->RC (-);[5]-->multiple mapping adjusted read count(-) <br>
	 */
	public static Map<String,double[]> getCountsBowtieOut(String bowtieFile){
		

		Map<String,double[]> back = new Hashtable<String,double[]>();
		Set<String> detected = new HashSet<String>();
		try {
			BufferedReader reader = new BufferedReader( new FileReader(bowtieFile));
			String line = null;

			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				double multmap = Double.parseDouble(f[6]);
				String[] f1 = f[0].split("#");
				double rc = 1d;
				if(f1.length >= 2){
					rc = Double.parseDouble(f1[1]);
				}
	
				String key = f[2];
				if(f.length > 10){
					key = f[9];
					if(detected.contains(f[0]+f[9])){
						continue;
					}
					else{
						detected.add(f[0]+f[9]);
					}
					
				}
				
				if(back.containsKey(key)){
					addToCount(f, back.get(key), rc, multmap);	
				}
				else{
					double[] t = new double[6];
					addToCount(f, t, rc, multmap);
					back.put(key, t);
				}


			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return back;
		
	}

	
    /**
     * This function makes non-redundant files for "tRNA" or "miR" (the different modes)
     * @param readAnnotationFile
     * @param out
     * @param group
     * @param orientation
     * @param totalMapped
     * @param mode
     */
    public static void getTRNAstat(Map<Integer,ReadData> readAnnotation, String out, String group, String orientation, double totalMapped, String mode){

    	Map<String,double[]> count = new Hashtable<String,double[]>();
    	try {


    		double libMapped = 0d;

    		for(int id : readAnnotation.keySet()){

 
    			for(AnnotData ad : readAnnotation.get(id).list){
    				if(ad.group.equals(group) && ad.orientation.equals(orientation)){

    					String key ="";
    					if(mode.equals("tRNA"))
    						key = ad.name.split("-")[1].split("=")[0].split("_")[0].split(":")[0];
    					else if(mode.equals("miR")){
    						key = ad.name.substring(4, ad.name.length());
    					}
    					else if(mode.equals("miRspec")){
    						key = ad.name;
    					}


    					libMapped += (double)readAnnotation.get(id).rc/(double)readAnnotation.get(id).list.size();

    					if(count.containsKey(key)){
    						count.get(key)[0] +=  1d/(double)readAnnotation.get(id).list.size();
    						count.get(key)[1] +=  (double)readAnnotation.get(id).rc/(double)readAnnotation.get(id).list.size();
    					}
    					else{
    						double[] tem = new double[2];
    						tem[0] = 1d/(double)readAnnotation.get(id).list.size();
    						tem[1] =(double)readAnnotation.get(id).rc/(double)readAnnotation.get(id).list.size();
    						count.put(key, tem);
    					}

    				}
    			}

    		}


    		BufferedWriter writer = new BufferedWriter (new FileWriter (out));
    		writer.write("antiCodon\tUR\tRC\tRC (adjusted)\tRPMlib\tRPMall\n");
    		for(String key : count.keySet()){
    			double rpm = count.get(key)[1]*1000000d/(double)libMapped;
    			double inputRPM = count.get(key)[1]*1000000d/(double)totalMapped;
    			writer.write(key+"\t"+count.get(key)[0]+"\t"+count.get(key)[1]+"\t"+"---"+"\t"+rpm+"\t"+inputRPM+"\n");
    		}

    		writer.close();
    		Sort.sortListBigger2Smaller(out, 2, true);
    	} catch (Exception e) {
    		// TODO Auto-generated catch block
    		e.printStackTrace();
    	}		
    }
    
	
	/**
	 * 
	 * @param readMap
	 * @param mode: regulates how the assignment of the reads will be performed: <br>
	 * mode=equal (read is assigned full to both strand if the reads maps to + and '-') <br>
	 * mode=sensePref (read is assigned to the + strand if it maps to both strands
	 * @return
	 */
	public static int[] getCountReadMap(Map<Integer,ReadData> readMap, String mode, String group){
		
		int[] back = new int[9];
		for(Integer readID : readMap.keySet()){

			int rc = readMap.get(readID).rc;
			boolean sense = false;
			boolean asense = false;
			for(AnnotData ad : readMap.get(readID).list){
				
				if(!(ad.group.equals(group)))
						continue;
				
				if(ad.orientation.equals("sense")){
					sense = true;
				}
				else if(ad.orientation.equals("antisense")){
					asense = true;
				}
				else{
					IO.warning("Did find unexpected strand in getCountReadMap(). libs.Stat");
				}
			}

			// assign fully to both strands if read maps to both strands!
			if(mode.equals("equal")){
				if(sense && asense){
					back[0]++; back[1]+=rc;back[2]++;back[3]+=rc;back[4]++;back[5]+=rc;
					back[6]+=(double)rc/(double)readMap.get(readID).multiple;
					back[7]+=(double)rc/(double)readMap.get(readID).multiple;
					back[8]+=(double)rc/(double)readMap.get(readID).multiple;
				}
				else if(sense){
					back[0]++; back[1]+=rc;back[2]++;back[3]+=rc;
					back[6]+=(double)rc/(double)readMap.get(readID).multiple;
					back[7]+=(double)rc/(double)readMap.get(readID).multiple;
				}
				else if(asense){
					back[0]++; back[1]+=rc;back[4]++;back[5]+=rc;
					back[6]+=(double)rc/(double)readMap.get(readID).multiple;
					back[8]+=(double)rc/(double)readMap.get(readID).multiple;
				}
			}
			// assign to sense strand
			else if(mode.equals("sensePref")){
				if(sense){
					back[0]++; back[1]+=rc;back[2]++;back[3]+=rc;
					back[6]+=(double)rc/(double)readMap.get(readID).multiple;
					back[7]+=(double)rc/(double)readMap.get(readID).multiple;

				}
				else if(asense){
					back[0]++; back[1]+=rc;back[4]++;back[5]+=rc;
					back[6]+=(double)rc/(double)readMap.get(readID).multiple;
					back[8]+=(double)rc/(double)readMap.get(readID).multiple;
				}
			}
			else{
				IO.warning("mode "+mode+" not found in getCountReadMap (libs.Stat)");
				System.exit(1);
			}
		}
		return back;
	}
	

		/**
		 * 
	 * takes a manipulated bowtie output file (*.gcross) and gives back the total read count for total / sense /antisense as a AlignData object
	 * @param file: can be either a *.parsed file or a *.gcross file
	 * @return
	 */
	public static AlignData profileFromBowtieOut(String file){

		double count[] = new double[9];
		Set<String> detected = new HashSet<String>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\\t");
				
				// the key: reference name + strand
				// each read is counted at most once for each strand
				String key = f[0]+f[1];
				double multmap = Double.parseDouble(f[6]);
				String[] f1 = f[0].split("#");
				double rc = 1d;
				if(f1.length >= 2){
					rc = Double.parseDouble(f1[1]);
				}
				/////////////////////////////////////////////
				/////////////////
				// for total count - count each read only once
				// check first if read was counted already
				if(detected.contains(f[0])){

				}
				else{
					detected.add(f[0]);
					count[6] ++;
					count[7] += rc;

				}			
				
				count[8] += rc/multmap;

				// if the read was counted before --> add new read count only
				// to the multiple mapping adjusted read count
				if(detected.contains(key)){

					addMultMappingAdjusted(f, count, rc, multmap);		

				}
				else{
					detected.add(key);
					addToCount(f, count, rc, multmap);

				}
			}
			reader.close();
		}catch (FileNotFoundException e){
			return null;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			return null;
		}
		AlignData ad = new AlignData(file,(int)count[6],(int)count[7],count[8],(int)count[0],(int)count[1],count[2],(int)count[3],(int)count[4],count[5]);
		return ad;		
	}
	
	
	/**
	 * add only the multiple adjusted to the count
	 * @param f
	 * @param count
	 * @param rc
	 * @param multmap
	 */
	private static void addMultMappingAdjusted(String[] f, double[] count, double rc, double multmap){
		

		// gcross format
		if(f.length >= 11){
			if(f[1].equals(f[10])){
				count[2] += rc/multmap;
			}
			else{
				count[5] += rc/multmap;
			}
		}
		// parsed format
		else{
			if(f[1].equals("+")){
				count[2] += rc/multmap;
			}
			else{
				count[5] += rc/multmap;
			}
		}
	}
	
	/**
	 * 
	 * @param f
	 * @param count
	 * @param rc
	 * @param multmap
	 */
	private static void addToCount(String[] f, double[] count, double rc, double multmap){

		// gcross format
		if(f.length >= 11){
			if(f[1].equals(f[10])){
				addToCount(count,0,rc,multmap);
			}
			else{
				addToCount(count,3,rc,multmap);
			}
		}
		// parsed format
		else{
			if(f[1].equals("+")){
				addToCount(count,0,rc,multmap);
			}
			else{
				addToCount(count,3,rc,multmap);
			}
		}
	}
	/**
	 * adds the read counts to the count array
	 * @param count
	 * @param index
	 * @param rc
	 * @param multmap
	 */
	private static void addToCount(double[] count, int index, double rc, double multmap){
		
		count[index]++;	
		count[index+1] += rc;
		count[index+2] += rc/multmap;
	}
	

	

	public static void  getProcessingStatSimple(List<MapData> mapDataList, BedDataRegion rd, int windowSize, Map<String,int[]> back){
		
//		Map<String,int[]> back = new Hashtable<String,int[]>();
		
		try {

			
			int regionLength = rd.end - rd.start + 1;
//			System.out.println("-------------"+rd.start+"-"+rd.end);
			for(MapData md : mapDataList){
			
//				System.out.println(md.start+"-"+md.end);
				
				int start = md.start; // difference of read to first position
				int end = regionLength - md.end;


				
				///////////////////////////////////////////////////////
				///
				
				if(start <= windowSize){
					Util.addIntMap(back, rd.name, md.count, 3, 0);
				}
				else if(end <= windowSize){
					Util.addIntMap(back, rd.name, md.count, 3, 1);
				}
				else{
					Util.addIntMap(back, rd.name, md.count, 3, 2);
				}


			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
//		return back;
	}
	
	

	
	public static Map<String,int[]> getProcessingStatSimple(String parsed, Map<String,Integer> lengthMap, int windowSize){
		
		Map<String,int[]> back = new Hashtable<String,int[]>();
		
		try {
			BufferedReader reader = new BufferedReader (new FileReader(parsed));

			String line = null;
			while((line = reader.readLine()) != null){
			
				String[] f = line.split("\\t"); // the line
				int rc = Integer.parseInt(f[0].split("#")[1]);
				
				String name = f[2];
				
				int start = 100; // difference of read to first position
				int end = 100; // difference to the 3p end
				// gcross format
				boolean full = true;
				///////////////////////////////////////////////////
				/////// gcross format
				//////
				if(f.length >= 10){
					name = f[9];
					// 	same strand
					if(f[1].equals(f[10])){
						if(f[1].equals("+")){
							start =   Integer.parseInt(f[3]) + 1 -  Integer.parseInt(f[11]) ; // number of bases from the beginning // 0-based
							end = Integer.parseInt(f[12]) - (Integer.parseInt(f[3]) + f[4].length()) ;  
						}
						else{
							start = Integer.parseInt(f[12]) - (Integer.parseInt(f[3]) + f[4].length());
							end = Integer.parseInt(f[3]) + 1 -  Integer.parseInt(f[11]);
//							System.out.println(start);
						}
					}
					// if read and tRNA are on different strands
					else{
						continue;
					}
				}
				
				///////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				// parsed format
				else{
					full = false;
					if(f[1].equals("+")){
						start = Integer.parseInt(f[3]);
						int len = 1000;
						if(lengthMap.containsKey(f[2])){
							len = lengthMap.get(f[2]);
						}
						else{
							IO.warning(f[2]+" not found in length map");
						}
//						System.out.println("length "+len);
						end = len - (start + f[4].length()) ;
					}
				}
				
				///////////////////////////////////////////////////////
				///
				
				if(start <= windowSize){
					Util.addIntMap(back, name, rc, 3, 0);
				}
				else if(end <= windowSize){
					Util.addIntMap(back, name, rc, 3, 1);
				}
				else{
					Util.addIntMap(back, name, rc, 3, 2);
				}


			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
		return back;
	}
	
	

	/**
	 * make the RNA composition as a function of read length
	 * @param readAnnot
	 * @param output
	 * @param referenceExpr
	 * @param readFileCol
	 * @param index
	 */
	public static void makeRNAcompositionReadLength(Map<Integer,ReadData> readAnnot, String output, 
			String readFile, int readFileCol, int index, boolean sensePreference){
		
		Map<Integer,Map<String,double[]>>  map = makeRNAcompositionLength(readAnnot, sensePreference);

		Map<Integer,Double> readLen = getReadLengthfromFile(readFile, readFileCol);
		makeRNAcompositionMatrixAsFunctionOfLength(map, output, index, readLen);
	}
	

	/**
	 * make the read length files of the different RNA types in the analysis 
	 * @param readAnnot
	 * @param output
	 * @param referenceExpr
	 * @param readFileCol
	 * @param index
	 */
	public static void makeRNAassignedReadLength(Map<Integer,ReadData> readAnnot, String outfolder, 
			int minReadLength, int maxReadLength, boolean sensePreference){
		
		new File(outfolder).mkdir();
		Map<Integer,Map<String,double[]>>  map = makeRNAcompositionLength(readAnnot,sensePreference);
		
		Map<String,double[][]> convMap = convertMap(map, maxReadLength);
		writeRNAreadLengths(convMap,outfolder, minReadLength, maxReadLength);

	}
	
	public static void writeRNAreadLengths(Map<String,double[][]> convMap , String outfolder,int minReadLength, int maxReadLength){
		
		for(String name : convMap.keySet()){
			double[] total = getCounts(convMap.get(name));
			String outfile = outfolder + File.separator + name+".readLen";
			try {
				BufferedWriter writer = new BufferedWriter ( new FileWriter(outfile));
				writer.write("Read Length (nt)\tUR\tPercentage_UR\tRC\tPercentage_RC\tRPM\n");
				for(int i = minReadLength; i <= maxReadLength; i++){
					double rcperc = 100d*convMap.get(name)[i][1]/total[1];
					double rpmperc = 100000d*convMap.get(name)[i][1]/total[1];
					double urperc = 100d*convMap.get(name)[i][0]/total[0];
					writer.write(i+"\t"+convMap.get(name)[i][0]+"\t"+urperc+"\t"+convMap.get(name)[i][1]+"\t"+rcperc+"\t"+rpmperc+"\n");
				}
				
				writer.close();
				if(GVars.graphics){
					graphics.getGraph(outfile, GVars.graphs+File.separator+Util.getFileBaseName(outfile)+".png", 1);
				}

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
	}
	
	private static double[] getCounts(double[][] countsReadLen){
		double[] back = new double[2];
		for(int i = 0; i < countsReadLen.length; i++){
			back[0]+= countsReadLen[i][0];
			back[1]+= countsReadLen[i][1];
		}
		return back;
	}
	
	private static Map<String,double[][]>  convertMap(Map<Integer,Map<String,double[]>>  map, int maxReadLength ){
		
		Map<String,double[][]> back = new Hashtable<String,double[][]>();
		for(Integer len : map.keySet()){
			Map<String,double[]> tempmap = map.get(len);
			/////////////////////////////////////////////////////////////
			// apply the classification as defined by the libsStrings 
			tempmap = applyLibStrings(tempmap);
			for(String key : tempmap.keySet()){
				if(back.containsKey(key)){
					back.get(key)[len][0] += tempmap.get(key)[0];
					back.get(key)[len][1] += tempmap.get(key)[1];
				}
				else{
					double[][] temp = new double[maxReadLength+1][2];
					temp[len][0] += tempmap.get(key)[0];
					temp[len][1] += tempmap.get(key)[1];
					back.put(new String(key), temp);
				}
			}
		}
		return back;
	}
	 
	/**
	 * read the information from a read length file
	 * @param readFile
	 * @param col
	 * @return
	 */
	public static Map<Integer,Double> getReadLengthfromFile(String readFile, int col){
		
		Map<Integer,Double> back = new Hashtable<Integer,Double>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(readFile));
			String line = reader.readLine();
			while((line = reader.readLine()) != null){
				String[] f = line.split("\\t");
				back.put(Integer.parseInt(f[0]), Double.parseDouble(f[col]));
			}
			
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return back;
	}

	/**
	 * 
	 * @param map
	 * @param output
	 * @param referenceExpr
	 * @param index
	 */
	public static void makeRNAcompositionMatrixAsFunctionOfLength(Map<Integer,Map<String,double[]>>  map, String output, int index, Map<Integer,Double> readLen ){
		
		try {
			Locale locale  = new Locale("en", "UK");
			String pattern = "###.##";

			DecimalFormat df = (DecimalFormat)
			        NumberFormat.getNumberInstance(locale);
			df.applyPattern(pattern);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(output));
			List<String> names = getNames(map);
			names.add("otherRNAs");
			writer.write("read length (nt)\tun-assigned");
			for(String name : names){
				writer.write("\t"+name);
			}
			writer.write("\n");
			for(Integer len : map.keySet()){
				
				double totalRC = 0;
				if(readLen.containsKey(len)){
					totalRC = readLen.get(len);
				}
				else{
					continue;
				}
				Map<String,double[]> countMap = map.get(len);
				
				double beforemapped = Util.getCountMapSum(countMap);
				if(index == 0)
				 beforemapped = Util.getCountMapSumUR(countMap);
				
				/////////////////////////////////////////////////////////////
				// apply the classification as defined by the libsStrings 
				countMap = applyLibStrings(countMap);

				
				double unAssigned = 100d*(totalRC - beforemapped)/totalRC;

				writer.write(len+"\t"+df.format(unAssigned));
				for(String name : names){
					double rc = 0;
					if(countMap.containsKey(name)){
						rc = countMap.get(name)[1];
						if(index == 0)
							rc = countMap.get(name)[0];
					}

					
					double perc = 100d*rc/totalRC;
					writer.write("\t"+df.format(perc));
				}
				writer.write("\n");
				
			}
			writer.close();
			Sort.sortListBigger2Smaller(output, 0, true);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	
	public static List<String> getNames(Map<Integer,Map<String,double[]>>  map){
		
		if(GVars.libsStringNames != null){
			String[] names = GVars.libsStringNames.split("\\|");
			List<String> back = new ArrayList<String>();
			for(String n : names){
				back.add(new String(n));
			}
			return back;
		}
		else{
			Set<String> t = new HashSet<String>();
			for(Integer I : map.keySet()){
				for(String key : map.get(I).keySet()){
					t.add(new String (key));
				}
			}
			List<String> back = new ArrayList<String>(t);
			return back;
		}
	}
	
	/**
	 * this function classifies the mappings like defined in the libsStringTypes and libsStringNames
	 * @param countMap
	 * @return
	 */
	public static Map<String,double[]> applyLibStrings(Map<String,double[]> countMap){
		
		Map<String,double[]> back = new Hashtable<String,double[]>();
		if(GVars.libsStringTypes != null && GVars.libsStringNames != null){
			String[] types = GVars.libsStringTypes.split("\\|");
			String[] names = GVars.libsStringNames.split("\\|");
		
			for(int i = 0; i < types.length; i++){
				// the different categories
				String[] cats = types[i].split(";");
				for(String cat : cats){
					if(countMap.containsKey(cat)){
						Util.addRCDouble(back,new String (names[i]), countMap.get(cat)[0], countMap.get(cat)[1]);
					}
				}
			}
			// all assigned reds
			double origRC = Util.getCountMapSum(countMap);
			double origUR = Util.getCountMapSumUR(countMap);
			// reads assigned to user defined categories
			double catRC = Util.getCountMapSum(back);
			double catUR = Util.getCountMapSumUR(back);
			
			double[] t = new double[2];
			t[0] = origUR - catUR;
			t[1] = origRC - catRC;
			back.put("otherRNAs", t);
			return back;
		}
		else{
			return countMap;
		}
	}
	
	public static int[] makeRNAcomposition(Map<Integer,ReadData> readAnnot, String output, int referenceUR, int referenceRC,
			boolean sensePreference){
		
		Map<String,double[]> map = makeRNAcomposition(readAnnot, sensePreference);
		map = applyLibStrings(map);
		int[] back = Write.writeCountMap(output, map, referenceUR,  referenceRC , "un-assigned");
		return back;
		
	}
	
	private static Map<String,double[]>  makeRNAcomposition(Map<Integer,ReadData> readAnnot, 
			boolean sensePreference){
		
		Map<String,double[]> map = new Hashtable<String,double[]>();
		for(Integer I : readAnnot.keySet()){
			
			ReadData rd = readAnnot.get(I);
			int numberOfSense = rd.numberOfSense();

			if(sensePreference && numberOfSense > 0){

					double rc = (double) readAnnot.get(I).rc / (double)numberOfSense;
					double ur = 1d / (double)numberOfSense;

					for(AnnotData ad : readAnnot.get(I).list){
						if(ad.orientation.equals("sense")){
							String[] f = ad.name.split(":");
							if(f.length > 1){
								Util.addRCDouble(map,f[1]+"#"+ad.orientation, ur, rc);
							}
							else{
								Util.addRCDouble(map,  ad.group+"#"+ad.orientation, ur, rc);
							}
						}
					}

			}
			else {
				double rc = (double) readAnnot.get(I).rc / (double)readAnnot.get(I).list.size();
				double ur = 1d / (double)readAnnot.get(I).list.size();
				
				for(AnnotData ad : readAnnot.get(I).list){
					String[] f = ad.name.split(":");
					if(f.length > 1){
						Util.addRCDouble(map,f[1]+"#"+ad.orientation, ur, rc);
					}
					else{
						Util.addRCDouble(map,  ad.group+"#"+ad.orientation, ur, rc);
					}
				}

			}
			
		}
		return map;
	}
	
	
	private static Map<Integer,Map<String,double[]>>   makeRNAcompositionLength(Map<Integer,ReadData> readAnnot,
			boolean sensePreference){

		//		Map<String,double[]> map = new Hashtable<String,double[]>();
		Map<Integer,Map<String,double[]>> map = new Hashtable<Integer,Map<String,double[]>>();
		for(Integer I : readAnnot.keySet()){

			ReadData rd = readAnnot.get(I);
			// the number of sense mappings
			int readLen = readAnnot.get(I).read.length();



			if(map.containsKey(readLen)){
				addAnnot(rd, map.get(readLen), readAnnot.get(I).rc, 
						sensePreference);
			}
			else{
				Map<String,double[]> t = new Hashtable<String,double[]>();
				addAnnot(rd, t, readAnnot.get(I).rc,sensePreference);
				map.put(readLen, t);
			}

		}


		return map;
	}
	
//	private static void addAnnot(List<AnnotData> list, Map<String,double[]> map, double rc){
//		
//		
//		
//		rc = rc / (double)list.size();
//		double ur = 1d / (double)list.size();
//		for(AnnotData ad : list){
//			String[] f = ad.name.split(":");
//			if(f.length > 1){
//				Util.addRCDouble(map,f[1]+"#"+ad.orientation, ur, rc);
//			}
//			else{
//				Util.addRCDouble(map,  ad.group+"#"+ad.orientation, ur, rc);
//			}
//		}
//	}

	
	private static void addAnnot(ReadData rd, Map<String,double[]> map, double rc, 
			boolean sensePreference){
		
		// number of sense mappings
		int numberOfSense = rd.numberOfSense();

		if(sensePreference && numberOfSense > 0){
			rc = rc / (double)numberOfSense;
			double ur = 1d / (double)numberOfSense;

			for(AnnotData ad : rd.list){
				if(ad.orientation.equals("sense")){
					
					String[] f = ad.name.split(":");
					if(f.length > 1){
						Util.addRCDouble(map,f[1]+"#"+ad.orientation, ur, rc);
					}
					else{
						Util.addRCDouble(map,  ad.group+"#"+ad.orientation, ur, rc);
					}
				}
			}
		}
		else{
		
			rc = rc / (double)rd.list.size();
			double ur = 1d / (double)rd.list.size();
			for(AnnotData ad : rd.list){
				String[] f = ad.name.split(":");
				if(f.length > 1){
					Util.addRCDouble(map,f[1]+"#"+ad.orientation, ur, rc);
				}
				else{
					Util.addRCDouble(map,  ad.group+"#"+ad.orientation, ur, rc);
				}
			}
		}
		


	}

	
	/**
	 * 
	 * @param groupedFile
	 * @param readAnnotFile
	 * @param output
	 * @param name
	 * @param orient
	 */
    public static void makeNonRedundant(String groupedFile,  Map<Integer,ReadData> readAnnotation, 
    		String output, String orient,
    		int libMapped, int totalMapped){


    	IO.log(GVars.logFile, 1, "Write out the non-redundant file of: "+output+" for orientation: "+orient, true);
		Map<String, Set<Integer>> map = ReadData.convertReadAnnotationShort(readAnnotation);
//		Map<String, Set<Integer>> map = ReadData.convertReadAnnotation(readAnnotation);
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(groupedFile));
			BufferedWriter writer = new BufferedWriter ( new FileWriter(output));
			reader.readLine();
			writer.write(Util.getGroupedSAHeader()+"\n");
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\\t");


				
				
				String key = f[0]+"#"+f[6].split(",")[0]+"#"+orient;
//				String key = f[0]+"#"+f[6]+"#"+orient;

				if(map.containsKey(key)){
						
					int rc = 0;
					int ur = 0;
					for(int index : map.get(key)){ // go over all indexes

						if(!(readAnnotation.get(index).asigned)){

							rc += readAnnotation.get(index).rc;
							ur++;
							readAnnotation.get(index).asigned = true;
						}
	
					
					}

					if(rc > 0){
						double librpm = 1000000d*rc/(double)libMapped;
						double totrpm = 1000000d*rc/(double)totalMapped;
						if(f.length >= 7)
							writer.write(f[0]+"\t"+ur+"\t"+rc+"\t"+f[2]+"\t"+librpm+"\t"+totrpm+"\t"+f[6]+"\n");
						else
							writer.write(f[0]+"\t"+ur+"\t"+rc+"\t"+f[2]+"\t"+librpm+"\t"+totrpm+"\n");
					}
					
				}
				else{
//					Helper.warning("Unexpected ");
				}
				

			}

			writer.close();
			reader.close();
			Sort.sortListBigger2Smaller(output, 2, true);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile,4,groupedFile+" or "+output+" not found",true);
			//			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block

			//			e.printStackTrace();
		}
	}



    public static void microRNAstat(String nonRedundantGroupedFile, String[] microRNAs, String top10Out, String speciesDistr, int nrTop){

    	if(new File(nonRedundantGroupedFile).exists()){
    		String tmp = nonRedundantGroupedFile+"tmp";
    		IO.copy(nonRedundantGroupedFile, tmp, false);
    		Sort.sortListBigger2Smaller(tmp, 2, true);
    		List<String> top10 = Read.readFileList(tmp, false);
    		if(top10.size() >= 1){
    			Write.writeString(top10Out, top10.get(0), false);
    			for(int i = 1; i < nrTop; i++){
    				if(top10.size() > i){
    					Write.writeString(top10Out, top10.get(i), true);
    				}
    			}
    		}
        	new File( tmp).delete();
    	}
    	// split microRNAs into species   	
    	if(microRNAs.length >=2){
    		Map<String,int[]> counts = countMicroRNAspecies(nonRedundantGroupedFile, microRNAs);
    		int miRmapped = 0;
    		for(String miR : counts.keySet())
    			miRmapped += counts.get(miR)[0];
    		
    		writeCounts(counts, speciesDistr, miRmapped);
    	}
    	

    }
    
    public static void writeCounts(Map<String,int[]> counts, String out, int microRNAmapped){
    	
    	Write.writeString(out, "species\tRC\tPercentage\tRPM", false);
		Locale locale  = new Locale("en", "UK");
		String pattern = "###.##";

		DecimalFormat df = (DecimalFormat)
		        NumberFormat.getNumberInstance(locale);
		df.applyPattern(pattern);
    	for(String miR : counts.keySet()){
    		double perc = 100d*(double)counts.get(miR)[0] / (double)microRNAmapped;
    		double rpm = 1000000d*(double)counts.get(miR)[0] / (double)microRNAmapped;
    		Write.writeString(out, miR+"\t"+counts.get(miR)[0]+"\t"+df.format(perc) +"\t"+df.format(rpm) , true);
    	}
    }
    
    public static Map<String,int[]> countMicroRNAspecies(String nonRedundantGroupedFile, String[] microRNAs){
    
    	
    	Map<String,int[]> back = new Hashtable<String,int[]>();
    	for(String miR : microRNAs)
    		back.put(miR, new int[1]);
    	
    	try {
			BufferedReader reader = new BufferedReader ( new FileReader(nonRedundantGroupedFile));
			String line = null;
			String header = reader.readLine();
			while(( line = reader.readLine()) != null){
				String[] f = line.split("\\t");
				
				for(String miR : microRNAs){
					if(line.startsWith(miR)){
						back.get(miR)[0] += Integer.parseInt(f[2]);
						break;
					}
				}
				
			}
			
			
			reader.close();
			return back;
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			IO.warning(nonRedundantGroupedFile+" not found");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	return back;
    }
	


}
