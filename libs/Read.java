package libs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import sRNAfuncTerms.Vars;



/**
 * functions to read file content
 * @author Michael
 *
 */
public class Read {

	/**
	 * Read ONE column of the file and give back a list 
	 * @param file: the file
	 * @param header: the file has a header?
	 * @param col  --> 0-based!!!
	 * @return
	 */
	public static Set<String> readFileColumn(String file, boolean header, int col){


		Set<String> back = new HashSet<String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header)
				reader.readLine();
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				if(col >= f.length ){
					IO.log(Vars.log, 4, "The column read from file: "+new File(file).getName()+" does not exist. Will quit now. ", true);
					IO.warning("The column read from file: "+new File(file).getName()+" does not exist. Will quit now. ");
					System.exit(1);
				}
				back.add(f[col].trim());
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return back;
	}


	/**
	 * gives back a map with Sequence <--> ID
	 * @param file
	 * @return
	 */
	public static Map<String,String> getFastaMapInverseSeed(String file, int seed){
		
		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null  ){
						String seq = sb.toString().toUpperCase().replace("U", "T").substring(0, seed);
						if(back.containsKey(seq)){
							String idt = back.get(seq) + "="+id.split("\\s+")[0];
							back.put(seq, idt);
						}
						else{
							back.put(seq,id.split("\\s+")[0]);					
						}
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();
				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null  ){
				String seq = sb.toString().toUpperCase().replace("U", "T").substring(0, seed);
				if(back.containsKey(seq)){
					String idt = back.get(seq) + "="+id.split("\\s+")[0];
					back.put(seq, idt);
				}
				else{
					back.put(seq,id.split("\\s+")[0]);					
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.noFile(GVars.logFile, file);
//			e.printStackTrace();
		}
		
		return back;
	}

	/**
	 * gives back a map with Sequence <--> ID
	 * @param file
	 * @return
	 */
	public static Map<String,String> getFastaMapInverseSeed(String file, int start, int end){
		
		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null  ){
						String seq = sb.toString().toUpperCase().replace("U", "T").substring(start-1, end);
						if(back.containsKey(seq)){
							String idt = back.get(seq) + "="+id.split("\\s+")[0];
							back.put(seq, idt);
						}
						else{
							back.put(seq,id.split("\\s+")[0]);					
						}
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();
				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null  ){
				String seq = sb.toString().toUpperCase().replace("U", "T").substring(start-1, end);
				if(back.containsKey(seq)){
					String idt = back.get(seq) + "="+id.split("\\s+")[0];
					back.put(seq, idt);
				}
				else{
					back.put(seq,id.split("\\s+")[0]);					
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.noFile(GVars.logFile, file);
//			e.printStackTrace();
		}
		
		return back;
	}

	
	public static String getHeader(String file){
		

		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			String line = reader.readLine();
			reader.close();
			return line;
		} catch (IOException e) {
			// TODO Auto-generated catch block
		}
	
		return null;
	}
	
	/**
	 * reads whole file content into a HAshSet
	 * the column index is 0-based!!!
	 * @param file
	 * @param keyCol: the column that should be used as key
	 * @param header
	 * @param the separater used in the file (\\t; \\s+ etc)
	 * @return
	 */
	public static Map<String,String> readFileMap(String file, int keyCol, boolean header, String sep){


		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header)
				reader.readLine();
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.replaceAll("\"", "").split(sep);
				back.put(f[keyCol], line.replaceAll("\"", "").trim());
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.writeToCommandLineL1(new File(file).getName()+" not found.");
			IO.log(GVars.logFile, 4, new File(file).getName()+" not found.", true);
			e.printStackTrace();
		}
		return back;
	}
	
	/**
	 * reads whole file content into a HAshSet
	 * the column index is 0-based!!!
	 * @param file
	 * @param keyCol: the column that should be used as key
	 * @param header
	 * @param the separater used in the file (\\t; \\s+ etc)
	 * @return
	 */
	public static Map<String,String> readFileMap(String file, int keyCol, int valueCol ,boolean header, String sep){


		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header)
				reader.readLine();
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.replaceAll("\"", "").split(sep);
				
				back.put(f[keyCol], f[valueCol]);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.writeToCommandLineL1(new File(file).getName()+" not found.");
			IO.log(GVars.logFile, 4, new File(file).getName()+" not found.", true);
			e.printStackTrace();
		}
		return back;
	}
	
	/**
	 * reads whole file content into a HAshSet
	 * the column index is 0-based!!!
	 * @param file
	 * @param keyCol: the column that should be used as key
	 * @param header
	 * @param the separater used in the file (\\t; \\s+ etc)
	 * @return
	 */
	public static Map<String,List<String>> readFileMapList(String file, int keyCol ,boolean header, String sep){


		Map<String,List<String>> back = new Hashtable<String,List<String>>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header)
				reader.readLine();
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.replaceAll("\"", "").split(sep);
				
				if(back.containsKey(f[keyCol])){
					back.get(f[keyCol]).add(line);
				}
				else{
					List<String> t = new ArrayList<String>(1);
					t.add(line);
					back.put(f[keyCol], t);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.writeToCommandLineL1(new File(file).getName()+" not found.");
			IO.log(GVars.logFile, 4, new File(file).getName()+" not found.", true);
			e.printStackTrace();
		}
		return back;
	}
	
	/**
	 * reads whole file content into a HAshSet
	 * the column index is 0-based!!!
	 * @param file
	 * @param keyCol: the column that should be used as key
	 * @param header
	 * @param the separater used in the file (\\t; \\s+ etc)
	 * @return
	 */
	public static Map<String,Set<String>> readFileMap_nonuniqueValues(String file, int keyCol, int valueCol ,boolean header, String sep){


		Map<String,Set<String>> back = new Hashtable<String,Set<String>>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header)
				reader.readLine();
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.replaceAll("\"", "").split(sep);
				if(back.containsKey(f[keyCol])){
					back.get(f[keyCol]).add(f[valueCol]);
				}
				else{
					Set<String> t = new HashSet<String>();
					t.add(f[valueCol]);
					back.put(f[keyCol], t);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.writeToCommandLineL1(new File(file).getName()+" not found.");
			IO.log(GVars.logFile, 4, new File(file).getName()+" not found.", true);
			e.printStackTrace();
		}
		return back;
	}

	
	/**
	 * reads whole file content into a List
	 * @param file
	 * @param header
	 * @return
	 */
	public static List<String> readFileList(String file, boolean header){


		List<String> back = new ArrayList<String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header){
				reader.readLine();
			}
			String line = null;
			while((line = reader.readLine()) != null){
				back.add(line.trim());
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.warning(file+" not found or not accessible");
			return null;
//			e.printStackTrace();
		}
		return back;
	}
	
	
	/**
	 * reads whole file content into a List
	 * @param file
	 * @param header
	 * @return
	 */
	public static List<String> readFileListColumn(String file, boolean header, int col, String sep){


		List<String> back = new ArrayList<String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader( file ));
			if(header)
				reader.readLine();
			String line = null;
			while((line = reader.readLine()) != null){
				String[] f = line.split(sep);
				back.add(f[col]);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.warning(file+" not found or not accessible");
			return null;
//			e.printStackTrace();
		}
		return back;
	}
	/**
	 * gives back a map with ID <--> Sequence length
	 * @param file
	 * @return
	 */
	public static Map<String,Integer> getFastaLengthMap(String file){
		
		Map<String,Integer> back = new Hashtable<String,Integer>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null  ){
						back.put(id.split("\\s+")[0], sb.length());
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();

				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null  ){
				back.put(id.split("\\s+")[0], sb.length());
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.warning(file+" not found in getFastaLengthMap");
//			e.printStackTrace();
		}
		
		return back;
	}
	
	/**
	 * gives back a map with Sequence <--> ID
	 * @param file
	 * @return
	 */
	public static Map<String,String> getFastaMapInverse(String file){
		
		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null  ){
						String seq = sb.toString().toUpperCase().replace("U", "T");
						if(back.containsKey(seq)){
							String idt = back.get(seq) + "="+id.split("\\s+")[0];
							back.put(seq, idt);
						}
						else{
							back.put(seq,id.split("\\s+")[0]);					
						}
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();
				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null  ){
				String seq = sb.toString().toUpperCase().replace("U", "T");
				if(back.containsKey(seq)){
					String idt = back.get(seq) + "="+id.split("\\s+")[0];
					back.put(seq, idt);
				}
				else{
					back.put(seq,id.split("\\s+")[0]);					
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Error while reading "+new File(file).getName()+". Will quit now. ", true);
			IO.warning("Error while reading "+new File(file).getName()+". Will quit now. ");
//			e.printStackTrace();
		}
		
		return back;
	}
	/**
	 * gives back a map with ID <--> Sequence. The sequence will be converted to DNA (in case it is RNA)
	 * @param file: a file in fasta format
	 * @return: 
	 */
	public static Map<String,String> getFastaMap(String file){
		
		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null  ){
						if(back.containsKey(id.split("\\s+")[0])){
							IO.warning(id+" is duplicated!!");
							IO.log(GVars.logFile, 3, id+" is duplicated and therefore ignored!", true);
						}
						else{
							back.put(new String(id.split("\\s+")[0]), sb.toString().toUpperCase().replace("U", "T"));
						}
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();

				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null  ){
				back.put(new String(id.split("\\s+")[0]), sb.toString().toUpperCase().replace("U", "T"));
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.warning(new File(file)+" not found (SeqUtil.getFastaMap). Will quit now!");
			IO.log(GVars.logFile, 4, new File(file)+" not found (SeqUtil.getFastaMap). Will quit now!", true);
			System.exit(1);
//			e.printStackTrace();
		}
		
		return back;
	}
	
	/**
	 * gives back a map with ID <--> Sequence. The sequence will be converted to DNA (in case it is RNA)
	 * @param file: a file in fasta format
	 * @return: 
	 */
	public static Map<String,BedDataRegion> getFastaMapRegion(String file){
		
		Map<String,BedDataRegion> back = new Hashtable<String,BedDataRegion>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null  ){
						if(back.containsKey(id.split("\\s+")[0])){
							IO.warning(id+" is duplicated in "+new File(file).getName());
							IO.log(GVars.logFile, 3, id+" is duplicated in "+new File(file).getName(), true);
						}
						else{
							BedDataRegion bd = new BedDataRegion(new String(id.split("\\s+")[0]), 1,  sb.toString().length(),new String(id.split("\\s+")[0])
							,0,"+",sb.toString().toUpperCase().replace("U", "T"));
							back.put(new String(id.split("\\s+")[0]), bd);
						}
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();

				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null  ){
				BedDataRegion bd = new BedDataRegion(new String(id.split("\\s+")[0]), 1,  sb.toString().length(),new String(id.split("\\s+")[0])
				,0,"+",sb.toString().toUpperCase().replace("U", "T"));
				back.put(new String(id.split("\\s+")[0]), bd);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
//			IO.warning(new File(file)+" not found (SeqUtil.getFastaMap). Will quit now!");
			IO.log(GVars.logFile, 4, new File(file)+" not found (SeqUtil.getFastaMap). Will quit now!", true);
			System.exit(1);
//			e.printStackTrace();
		}
		
		return back;
	}
	
	
	
	
	/**
	 * filters out the microRNAs that the user whats to analyse
	 * @param in
	 * @param out
	 * @param species
	 */
	public static int getSpeciesMicroRNAs(String in, String out, String[] species){
		
		int c = 0;
		try {
			BufferedReader reader = new BufferedReader (new FileReader(in));
			BufferedWriter writer = new BufferedWriter (new FileWriter(out));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null && Util.checkSpecies(id, species)){
						c++;
						writer.write(">"+id.split("\\s+")[0]+"\n");
						writer.write(sb.toString().toUpperCase().replace("U", "T")+"\n");
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();

				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null && Util.checkSpecies(id, species)){
				c++;
				writer.write(">"+id.split("\\s+")[0]+"\n");
				writer.write(sb.toString().toUpperCase().replace("U", "T")+"\n");
			}
			reader.close();
			writer.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			IO.noFile(GVars.logFile, in);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.noAccess(GVars.logFile, out);

		}
		return c;
	}

	
	/**
	 * filters out the microRNAs that the user whats to analyse
	 * @param in
	 * @param out
	 * @param species
	 */
	public static int getSpeciesMicroRNAs(String in, String out, String[] species, int maxLength){
		
		int c = 0;
		try {
			BufferedReader reader = new BufferedReader (new FileReader(in));
			BufferedWriter writer = new BufferedWriter (new FileWriter(out));
			
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					if(id != null && Util.checkSpecies(id, species)){
						if(sb.length() <= maxLength){
							c++;
							writer.write(">"+id.split("\\s+")[0]+"\n");
							writer.write(sb.toString().toUpperCase().replace("U", "T")+"\n");
						}
						else{
							IO.log(GVars.logFile, 3, "Eliminate a read longer than "+maxLength+" nt "+id.split("\\s+")[0], true);
//							IO.warning();
						}
					}
					id = line.trim().replace(">", "");
					sb = new StringBuilder();

				}
				else{
					sb.append(line.trim());
				}
			}
			
			if(id != null && Util.checkSpecies(id, species)){
				c++;
				writer.write(">"+id.split("\\s+")[0]+"\n");
				writer.write(sb.toString().toUpperCase().replace("U", "T")+"\n");
			}
			reader.close();
			writer.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			IO.noFile(GVars.logFile, in);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.noAccess(GVars.logFile, out);

		}
		return c;
	}

	
	/**
	 * Reads a file in JSON format. !!! the file cannot contain lists !!!!
	 * @param file
	 * @return
	 */
	public static Map<String,String> readJsonSimple(String file){
		
		Map<String,String> back = new Hashtable<String,String>();
		try {
			BufferedReader reader = new BufferedReader (new FileReader(file));
			String line = null;
			while((line = reader.readLine()) != null){
				if(line.contains(":")){
					String[] f = line.split(":");
					if(f.length == 2){
						String key = f[0].replace("\"", "").trim();
						String val = f[1].replace("\"", "").replace(",", "").trim();
						back.put(key, val);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}
		return back;
	}


}
