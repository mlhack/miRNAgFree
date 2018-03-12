package libs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import sequences.FastaByteUC;
import sequences.MFasta;



public class Util {

	
	
	/**
	 * gives back the root name (cut of the extension). this function is important for filenames with more than one '.' in the name
	 * @param filename
	 * @return
	 */
	public static String getFileBaseName(String filename){
		
		String fn = new File(filename).getName();
		String[] f = fn.split("\\.");
		StringBuilder sb = new StringBuilder();
		sb.append(f[0]);
		for(int i = 1; i < f.length-1; i++){
			sb.append("."+f[i]);
		}
		return sb.toString();
		
	}
	
	/**
	 * Add data to the parameter Map
	 * @param info
	 * @param tag
	 * @param value
	 */
	public static void setInfo(Map<String,List<String>> info,String tag, String value){

		if(info.containsKey(tag)){
			info.get(tag).add(value);
		}
		else{
			List<String> temp = new ArrayList<String>(3);
			temp.add(value);
			info.put(tag, temp);
		}

	}
	/**
	 * 
	 * @param set
	 * @return
	 */
	static public Set<Thread> checkTreads(Set<Thread> set){
		
		Set<Thread> back = new HashSet<Thread>(set.size());
		for(Thread thread : set){
			if(thread.isAlive()){
				back.add(thread);
			}
		}
		return back;
	}
	
	// check if the id corresponds to any of the species 
	public static boolean checkSpecies(String id, String[] species){
		
		for(String spec : species){
			if(id.startsWith(spec)){
				return true;
			}
			if(spec.equals("all")){
				return true;
			}
			
		}
		return false;
	}
	
	
	/**
	 * gets back a random string
	 * @param length
	 * @return
	 */
	public static String getRandomString(int length){
		
		char[] chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ".toCharArray();
		StringBuilder sb = new StringBuilder();
		Random random = new Random();
		for (int i = 0; i < length; i++) {
		    char c = chars[random.nextInt(chars.length)];
		    sb.append(c);
		}
		String output = sb.toString();
		return output;
	}
	
	/**
	 * 
	 * @return: the header of all grouped files
	 */
	public static String getGroupedHeader(){
		return "name\tunique reads\tread count\tread count (mult. map. adj.)\tRPM (lib)\tRPM (total)\tcoordinateString\tRPM_adj (lib)\tRPM_adj (total)";
	}

	public static String getGroupedSAHeader(){
		return "name\tunique reads\tread count (SA)\tread count (MA)\tRPM (lib)\tRPM (total)\tcoordinateString\tRPM_adj (lib)\tRPM_adj (total)";
	}
	
	
	/**
	 * removes certain reads from the reads.fa input file
	 * @param readFile
	 * @param set
	 */
	public static void removeReadsFasta(String readFile, Set<String> set){

		try {
			BufferedReader reader = new BufferedReader (new FileReader(readFile));
			BufferedWriter writer = new BufferedWriter (new FileWriter(readFile+"T"));

			String line = null;
			while((line = reader.readLine()) != null){
				if(line.startsWith(">")){

					String seq = reader.readLine();
					String id = line.replace(">", "");
					if(set.contains(id)){

					}
					else{
						writer.write(line+"\n");
						writer.write(seq+"\n");
					}
				}

			}
			reader.close();
			writer.close();
			new File(readFile+"T").renameTo(new File(readFile));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	/**
	 * removes certain reads from the reads.fa input file
	 * @param readFile
	 * @param set
	 */
	public static void removeReadsFasta(String readFile, Map<Integer,ReadData> readAnnotation){

		try {
			BufferedReader reader = new BufferedReader (new FileReader(readFile));
			BufferedWriter writer = new BufferedWriter (new FileWriter(readFile+"T"));

			String line = null;
			while((line = reader.readLine()) != null){
				if(line.startsWith(">")){

					String seq = reader.readLine();
					String id = line.replace(">", "");
					String[] f = id.split("\\s+")[0].split("#");
					if(readAnnotation.containsKey(Integer.parseInt(f[0]))){

					}
					else{
						writer.write(line+"\n");
						writer.write(seq+"\n");
					}
				}

			}
			reader.close();
			writer.close();
			new File(readFile+"T").renameTo(new File(readFile));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	/**
	 * Write the UTR file into smaller chunks
	 * @param file
	 * @param outBase
	 * @param chunkSize: number of sequences in one chunk
	 * @return
	 */
	public static List<String> makeChunks(String file, String outBase, int chunkSize){

		List<String> back = new ArrayList<String>();
		MFasta mf = new MFasta(file);
		int c = 1; // chunk number counter
		int s = 0; // chunk size counter
		Map<String,String> tmpMap = new Hashtable<String,String>();
		while(mf.hasNext()){
			FastaByteUC data = mf.next();
			if(tmpMap.containsKey(data.getId())){
				IO.log(GVars.logFile, 3, "Duplicated entries "+data.getId(), true);
			}
			else{
				try {
					tmpMap.put(data.getId(), data.getSequenceString());
				} catch (UnsupportedEncodingException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

		String out = outBase+File.separator+"tmp_"+c+".txt";
		new File(out).delete();
		for(String id : tmpMap.keySet()){
			if(tmpMap.get(id).length() >= 25)
				Write.writeFasta(out, tmpMap.get(id), id,"", true);
			s++;
			if(s >= chunkSize){
				s= 0;
				c++;
				back.add(out);
				out = outBase+File.separator+"tmp_"+c+".txt";
				new File(out).delete();
			}
		}
		if(new File(out).isFile())
			back.add(out);
		return back;
	}

	
	/**
	 * 
	 * @param file
	 */
	public static List<Sort> getSRNAbenchFormatSortList(String file){
		
		try {
			List<Sort> list = new ArrayList<Sort>();
			BufferedReader reader = new BufferedReader( new FileReader(file));
			String line = null;
			while((line = reader.readLine()) != null){
				String seq = reader.readLine();
				Sort s = new Sort(Double.parseDouble(line.split("#")[1]),seq);
				list.add(s);
			}
			reader.close();
			Sort.sortBiggerToSmaller(list);
			return list;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * sorts a internal sRNAbench fasta format file (reads.fa)<br>
	 * Format:<br>
	 * >1#4  --> (ID#Read Count)
	 * AGCGGCAATGTACCAT
	 * @param file
	 */
	public static void sortSRNAbenchFormat(String file){
		
		try {
			List<Sort> list = new ArrayList<Sort>();
			BufferedReader reader = new BufferedReader( new FileReader(file));
			String line = null;
			while((line = reader.readLine()) != null){
				String seq = reader.readLine();
				double sortV = 1d;
				String[] f = line.split("#");
				if(f.length > 1){
					sortV = Double.parseDouble(f[1]);
				}
				Sort s = new Sort(sortV,line+"\t"+seq);
				list.add(s);
			}
			reader.close();
			Sort.sortBiggerToSmaller(list);
			BufferedWriter writer = new BufferedWriter( new FileWriter(file));
			for(Sort s : list){
				String[] f = s.line.split("\t");
				writer.write(f[0]+"\n");
				writer.write(f[1]+"\n");
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public static String getString(String[] f, String sep){
		
		if(f.length > 0){
			StringBuilder sb = new StringBuilder();
			sb.append(f[0]);
			for(int i = 1; i < f.length; i++){
				sb.append(sep+f[i]);
			}
			return sb.toString();
		}
		return "";
	}
	
	/**
	 * gives back a String with 'nr' of times repeated 'fill'
	 * @param nr
	 * @return
	 */
	public static String getCharString(int nr, char fill){
		
		StringBuilder sb = new StringBuilder();
		for(int i = 1; i <= nr; i++){
			sb.append(fill);
		}
		return sb.toString();
	}

	/**
	 * gives back a String with 'nr' of times repeated 'fill'
	 * @param nr
	 * @return
	 */
	public static String getString(int nr, String fill){
		
		StringBuilder sb = new StringBuilder();
		for(int i = 1; i <= nr; i++){
			sb.append(fill);
		}
		return sb.toString();
	}

	/**
	 * gives back a String with 'nr' spaces
	 * @param nr
	 * @return
	 */
	public static String getVoidString(int nr){
		
		StringBuilder sb = new StringBuilder();
		for(int i = 1; i <= nr; i++){
			sb.append(" ");
		}
		return sb.toString();
	}



	/**
	 * Supposes a map like: String <-->int[]. sums to an existing key or put a new key
	 * @param map: 
	 * @param key
	 * @param value
	 */
	public static void addIntMap(Map<String,int[]> map, String key, int value){
		
		if(map.containsKey(key.trim())){
			map.get(key.trim())[0] += value;
		}
		else{
			int[] t = new int[1];
			t[0]=value;
			map.put(new String(key).trim(), t);

		}
	}
	
	public static void addIntMapRC(Map<String,int[]> map, String key, int value){
		
		if(map.containsKey(key.trim())){
			map.get(key.trim())[0] ++;
			map.get(key.trim())[1] += value;
		}
		else{
			int[] t = new int[2];

			t[0]=1;
			t[1]=value;
			map.put(new String(key).trim(), t);

		}
	}
	
	/**
	 * Supposes a map like: String <-->int[]. sums to an existing key or put a new key
	 * @param map: 
	 * @param key
	 * @param value
	 */
	public static void addIntMap(Map<String,int[]> map, String key, int value, int size, int index){
		
		if(map.containsKey(key)){
			map.get(key)[index] += value;
		}
		else{
			int[] t = new int[size];
			t[index]=value;
			map.put(new String(key), t);

		}
	}
	
	/**
	 * adds Read count to a map
	 * @param map
	 * @param key
	 * @param value
	 */
	public static void addRC(Map<String,double[]> map, String key, int value){
		if(map.containsKey(key)){
			map.get(key)[0]++;
			map.get(key)[1]+=value;
		}
		else{
			double[] t = new double[2];
			t[0]=1;
			t[1]=value;
			map.put(key, t);
		}
	}
	/**
	 * adds Read count to a map
	 * @param map
	 * @param key
	 * @param value
	 */
	public static void addRCDouble(Map<String,double[]> map, String key, double ur, double rc){
		if(map.containsKey(key)){
			map.get(key)[0]+=ur;
			map.get(key)[1]+=rc;
		}
		else{
			double[] t = new double[2];
			t[0]=ur;
			t[1]=rc;
			map.put(new String(key), t);
		}
	}
	
	/**
	 * get a string representation of the Set separating the entries by '-'
	 * @param set
	 * @return
	 */
	public static String getSetString(Set<String> set){
		
		StringBuilder sb = new StringBuilder();
		for(String s : set)
			sb.append(s+"-");
		
		return sb.toString().substring(0, sb.length()-1);
	}
	
	
	/**
	 * 
	 * @param map
	 * @return
	 */
	public static double getCountMapSum(Map<String,double[]> map){
		
		double back=0d;
		for(String key : map.keySet()){
			back += map.get(key)[1];
		}
		return back;
	}
	public static double getCountMapSumUR(Map<String,double[]> map){
		
		double back=0d;
		for(String key : map.keySet()){
			back += map.get(key)[0];
		}
		return back;
	}
	
	
	
	
	
	/**
	 * check if bowtie index exists
	 * @param baseName
	 * @return
	 */
	public static boolean checkIndex(String baseName){
	
		
		int count = getNumberOfIndexFiles(baseName, "ebwt");
		if(count == 6){
			return true;
		}
		else if(count == 0){
			int countB = getNumberOfIndexFiles(baseName, "ebwtl");
			if(countB == 6){
				return true;
			}
			else if(countB == 0){
				return false;
			}
			else{
				IO.warning("The bowtie index basename: "+baseName+" is "
						+ "incomplete (not all 6 files exist)! \n"
						+ " **** Please try running the populate tool again, or relaunch bowtie-build! **** "
						+ "\n Will quit now.");
				System.exit(0);
	
			}

		}
		else{
			IO.warning("The bowtie index basename: "+baseName+" is "
					+ "incomplete (not all 6 files exist)! \n"
					+ " **** Please try running the populate tool again, or relaunch bowtie-build! **** "
					+ "\n Will quit now.");
			System.exit(0);
		}

		return false;
	}
	
	/**
	 * 
	 * count the number of bowtie index files
	 * @param baseName
	 * @param indexType
	 * @return
	 */
	private static int getNumberOfIndexFiles(String baseName, String indexType){
		
		int count = 0;
		// check normal indexes
		if(new File(baseName+".1."+indexType).exists()){
			count++;
		}
		if(new File(baseName+".2."+indexType).exists()){
			count++;
		}
		if(new File(baseName+".3."+indexType).exists()){
			count++;
		}
		if(new File(baseName+".4."+indexType).exists()){
			count++;
		}
		if(new File(baseName+".rev.1."+indexType).exists()){
			count++;
		}
		if(new File(baseName+".rev.2."+indexType).exists()){
			count++;
		}
		
		return count;
	}
	
	
	/**
	 * 
	 * @param parameterString
	 * @param description
	 * @param maxParameterLen
	 * @param totalWindowlen
	 */
	public static void printArguments(String parameterString, String description, int maxParameterLen, int totalWindowlen){
		int realParameterLen = parameterString.length();
		int fillUp = maxParameterLen - realParameterLen;
		if(fillUp <= 0)
			fillUp = 1;
		
		String fillString = Util.getCharString(fillUp, ' ');
		
		System.out.print(parameterString+fillString);
		
		int descriptionLen = totalWindowlen - maxParameterLen;
		
		List<String> descChunks = splitString(description, descriptionLen);
		System.out.println(descChunks.get(0));//.replace("<>", "\n"));

		String secondFillString = Util.getCharString(maxParameterLen, ' ');

		for(int i = 1; i < descChunks.size(); i++){
			System.out.println(secondFillString+descChunks.get(i));//.replace("<>", "\n"));
		}
		
		
	}

	/**
	 * 
	 * @param str
	 * @param chunkLen
	 * @return
	 */
	private static List<String> splitString(String str, int chunkLen){
		
		List<String> back = new ArrayList<String>();
		
		String[] splitted = str.split("\\s+");
		StringBuilder sb = new StringBuilder();
		for(String strtemp : splitted){
			sb.append(strtemp+" ");
			if(sb.length() > chunkLen){
				back.add(sb.toString());
				sb = new StringBuilder();
			}
		}
		back.add(sb.toString());
		return back;
	}
	
	
	/**
	 * get a string representation of the Set separating the entries by '-'
	 * @param set
	 * @param the separator
	 * @return
	 */
	public static String getSetString(Set<String> set, String sep){
		
		StringBuilder sb = new StringBuilder();
		for(String s : set)
			sb.append(s+sep);
		
		return sb.toString().substring(0, sb.length()-sep.length());
	}
	/**
	 * get a string representation of the Set separating the entries by '-'
	 * @param set
	 * @param the separator
	 * @return
	 */
	public static String getListString(List<String> set, String sep){
		
		StringBuilder sb = new StringBuilder();
		for(String s : set)
			sb.append(s+sep);
		
		return sb.toString().substring(0, sb.length()-sep.length());
	}
}
