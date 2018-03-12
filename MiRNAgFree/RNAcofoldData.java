package miRNAgFreeGit;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

import libs.GVars;
import libs.IO;
import libs.Util;
import sequences.SeqUtil;


/**
 * class to store parsed RNAcofold output
 * @author MH
 *
 */
public class RNAcofoldData {

	public String sec5p;
	public String sec3p;
	String struc5p;
	String struc3p;
	double energy;
	double energyRatio;
	public String name;
	int bindings;
	int trailingNonBindingsTotal;
	GFcluster clust5p;
	GFcluster clust3p;
	int droshaDicerFluc;
	
	public RNAcofoldData(String secs, String structure, String name){
		
		String[] f = secs.split("&");
		this.sec5p = f[0].replaceAll("U", "T");
		this.sec3p = f[1].replaceAll("U", "T");
		String[] st1 = structure.split("\\s+");
//		System.out.println(structure);
		this.energy = -1d*Double.parseDouble(structure.split("-")[1].trim().replaceAll("\\)", ""));
		String[] st2 = st1[0].split("&");
		this.struc5p = st2[0];
		this.struc3p = st2[1];
		this.name = name.replaceAll(">", "");
		this.energyRatio = 2d*this.energy/(double)(this.sec5p.length() + this.sec3p.length());
		this.calcBindings();
		int[] nb3p =  this.getTrailing3pNonBindings();
		int[] nb5p = this.getTrailing5pNonBindings();
		this.trailingNonBindingsTotal = nb3p[0]+nb3p[1]+nb5p[0]+nb3p[1];
		this.droshaDicerFluc = getDroshaDicerFluctuation();
	}
	
	public RNAcofoldData(){
		
	}
	
	
	/**
	 * get the number of bindings
	 * 
	 */
	private void calcBindings(){
		int count = 0;
		for(char c : this.struc5p.toCharArray()){
			if(c == '(')
				count++;
		}
		
		this.bindings =  count;
	}
	
	/**
	 * get the number of bindings
	 * 
	 */
	public int getBindings(){
		return this.bindings;
	}
	
	/**
	 * [0] --> non-bindigs in 5p; [1] --> non-bindings in 3p
	 * @return
	 */
	public int[] getNonBinding(){
		
		int[] back = new int[2];
		back[0] = this.sec5p.length() - this.bindings;
		back[1] = this.sec3p.length() - this.bindings;
		return back;
	}
	
	/**
	 * the number of non-binding position at the 3' end for [0] 5p arm and [1] 3p arm
	 * @return
	 */
	public int[] getTrailing3pNonBindings(){
		
		int[] back = new int[2];
		char[] spl = this.struc5p.toCharArray();
		for(int i = spl.length-1; i >0 ; i--){
			if(spl[i] == '.'){
				back[0]++;
			}
			else
				break;
		}

		spl = this.struc3p.toCharArray();
		for(int i =spl.length -1; i >= 0; i--){
			if(spl[i] == '.'){
				back[1]++;
			}
			else
				break;
		}
		
		return back;
	}
	/**
	 * the number of non-binding position at the 5' end for [0] 5p arm and [1] 3p arm
	 * @return
	 */
	public int[] getTrailing5pNonBindings(){
		
		int[] back = new int[2];
		char[] spl = this.struc5p.toCharArray();
		for(int i = 0; i < spl.length ; i++){
			if(spl[i] == '.'){
				back[0]++;
			}
			else
				break;
		}

		spl = this.struc3p.toCharArray();
		for(int i = 0; i < spl.length ; i++){
			if(spl[i] == '.'){
				back[1]++;
			}
			else
				break;
		}
		
		return back;
	}
	
	public int[][] getInnerNonbindingsLength(){
		
		int[][] back = new int[2][4];
		String[] spl = this.struc5p.split("\\(");
//		System.out.println(this.struc5p);
		for(int i = 1; i < spl.length -1; i++){
//			System.out.println(spl[i]);
			if(spl[i].length() == 1){
				back[0][0]++;
			}
			if(spl[i].length() == 2){
				back[0][1]++;
			}
			if(spl[i].length() == 3){
				back[0][2]++;
			}
			if(spl[i].length() > 3){
				back[0][3]++;
			}

		}

		spl = this.struc3p.split("\\)");
		for(int i = 1; i < spl.length -1; i++){
			if(spl[i].length() == 1){
				back[1][0]++;
			}
			if(spl[i].length() == 2){
				back[1][1]++;
			}
			if(spl[i].length() == 3){
				back[1][2]++;
			}
			if(spl[i].length() > 3){
				back[1][3]++;
			}

		}

		return back;
	}
	
	public int[] getInnerNonbindings(){
		
		int[] back = new int[2];
		String[] spl = this.struc5p.split("\\(");
		for(int i = 1; i < spl.length -1; i++){
			if(spl[i].contains(".")){
				back[0]+=spl[i].length();
			}
		}

		spl = this.struc3p.split("\\)");
		for(int i = 1; i < spl.length -1; i++){
			if(spl[i].contains(".")){
				back[1]+=spl[i].length();
			}
		}

		return back;
	}
	
	/**
	 * checks if the hybrid has inner bulbs, i.e. if the hybrid is 'perfect'
	 * @return
	 */
	public boolean hasInnerNonBindings(){
		
		String[] spl = this.struc5p.split("\\(");
		for(int i = 1; i < spl.length -1; i++){
			if(spl[i].contains(".")){
				return true;
			}
		}

		spl = this.struc3p.split("\\)");
		for(int i = 1; i < spl.length -1; i++){
			if(spl[i].contains(".")){
				return true;
			}
		}
		return false;
	}
	///////////////////////////////////////////////////////////
	/// STATIC FUNCTIONS

	/**
	 * 
	 * @param struc
	 * @param split
	 * @return // 0--> total number of non-bindings ; 1 --> maximum length of bulge or internal loop; 2 --> number of bulges or loops
	 */
	public static int[] getStrucProperties(String struc,String split){
		
		int[] back = new int[3]; 

		String[] t = struc.split("\\"+split+"+");
//		System.out.println(struc);
		back[1] = -1;
		back[2] = t.length;
		for(String o : t){
			back[0]+=o.length();
			if(o.length() > back[1]){
				back[1] = o.length();
			}
		}
		return back;
	}

	/**
	 * 
	 * @param list
	 * @return list of RNAcofoldData objects 
	 */
	public static List<RNAcofoldData> getCofoldDataLax(List<String> list){
		
		List<RNAcofoldData> back = new ArrayList<RNAcofoldData>();
		
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");
//					System.out.println(list.get(i+2));
					int[] trailing1 = getTrailingNonBinding(f[0].toCharArray());
					int[] trailing2 = getTrailingNonBinding(f[1].toCharArray());
					int[] b5p = bindings(f[0]);
					int[] b3p = bindings(f[1]);
					if( ((b5p[0] == 0 || b5p[1] == 0) && (b3p[0] == 0 || b3p[1] == 0)) && trailing1[1] - trailing2[0] >= 1 && trailing1[1] - trailing2[0] <= 3 
							&& trailing2[1] - trailing1[0] >= 1 && trailing2[1] - trailing1[0] <= 3){
						RNAcofoldData data = new RNAcofoldData(new String(list.get(i+1)), new String (list.get(i+2)), new String (list.get(i)));
						back.add(data);
					}
				}
			}
		}
		
		return back;
	}
	
	/**
	 * gives back the total number of bases that deviate from a perfect 2nt 3' overhang pattern
	 * @return
	 */
	public int getDroshaDicerFluctuation(){
		
		int[] trailing1 = getTrailingNonBinding(this.struc5p.toCharArray());
		int[] trailing2 = getTrailingNonBinding(this.struc3p.toCharArray());
		
		return (   Math.abs(2- (trailing1[1] - trailing2[0])) + Math.abs (2 - (trailing2[1] - trailing1[0])) );
	}
	
	public static List<RNAcofoldData> getCofoldDataStrict(List<String> list){
		
		List<RNAcofoldData> back = new ArrayList<RNAcofoldData>();
		
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");
//					System.out.println(list.get(i+2));
					int[] trailing1 = getTrailingNonBinding(f[0].toCharArray());
					int[] trailing2 = getTrailingNonBinding(f[1].toCharArray());
					int[] b5p = bindings(f[0]);
					int[] b3p = bindings(f[1]);
					if(((b5p[0] == 0 || b5p[1] == 0) && (b3p[0] == 0 || b3p[1] == 0)) 
							&& trailing1[1] - trailing2[0] == 2 
							&& trailing2[1] - trailing1[0] == 2){
						RNAcofoldData data = new RNAcofoldData(new String(list.get(i+1)), new String (list.get(i+2)), new String (list.get(i)));
						back.add(data);
					}
				}
			}
		}
		return back;
	}
	

	public static List<RNAcofoldData> getCofoldDataMinBinding(List<String> list, int minBinding){
		
		List<RNAcofoldData> back = new ArrayList<RNAcofoldData>();
		
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");
					//					System.out.println(list.get(i+2));
					int[] b5p = bindings(f[0]);
					int[] b3p = bindings(f[1]);
					// 
					if( b5p[0] >= minBinding  && ((b5p[0] == 0 || b5p[1] == 0) && (b3p[0] == 0 || b3p[1] == 0))){

						RNAcofoldData data = new RNAcofoldData(new String(list.get(i+1)), new String (list.get(i+2)), new String (list.get(i)));
						back.add(data);
					}
				}
			}
		}

		return back;
	}

	/**
	 * 
	 * @param struc
	 * @return
	 */
	public static int[] bindings(String struc){
		
		int[] count = new int[2];
		for(char c : struc.toCharArray()){
			if(c == '(')
				count[0]++;
			else if(c==')'){
				count[1]++;
			}
		}
		return count;
	}
	
	/**
	 * get only Drosha/Dicer duplexes back
	 */
	public static List<String> parseStrict(List<String> list){

		List<String> back = new ArrayList<String>();
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");
//					System.out.println(list.get(i+2));
					int[] trailing1 = getTrailingNonBinding(f[0].toCharArray());
					int[] trailing2 = getTrailingNonBinding(f[1].toCharArray());
					if(trailing1[1] - trailing2[0] == 2 && trailing2[1] - trailing1[0] == 2){
//					if(f[0].startsWith("(") && f[1].startsWith(")") && f[0].endsWith("..") && f[1].endsWith("..") && !(f[0].endsWith("...")) && !(f[1].endsWith("..."))){
						back.add(list.get(i));
						back.add(list.get(i+1));
						back.add(list.get(i+2));
					}
				}
			}
		}

		return back;
	}
	
	/**
	 * takes a list of 
	 * @param list
	 * @return
	 */
	public static Map<String,int[]> parseToPattern(List<String> list,int casseteSize){
		
		Map<String,int[]> countMap = new Hashtable<String,int[]>();
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");
					f[0] = f[0].substring(casseteSize,f[0].length());
					f[0] = f[0].substring(0,f[0].length()-casseteSize);
					
					f[1] = f[1].substring(casseteSize,f[1].length());
					f[1] = f[1].substring(0,f[1].length()-casseteSize);
					
					add(countMap, "5p:"+f[0], 1);
					add(countMap, "3p:"+f[1], 1);
					add(countMap, "whole:"+f[0]+"&"+f[1], 1);
				}
			}
		}
		return countMap;
	}
	/**
	 * takes a list of 
	 * @param list
	 * @return
	 */
	public static Map<String,int[]> parseToPattern(List<String> list,Map<String,StringBuilder> microRNAmap){
		
		Map<String,int[]> countMap = new Hashtable<String,int[]>();
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");

					
					add(countMap, "5p:"+f[0], 1);
					add(countMap, "3p:"+f[1], 1);
					add(countMap, "whole:"+f[0]+"&"+f[1], 1);
					addMiR(microRNAmap,"5p:"+f[0],list.get(i).replace(">", ""));
					addMiR(microRNAmap,"3p:"+f[1],list.get(i).replace(">", ""));
					addMiR(microRNAmap,"whole:"+f[0]+"&"+f[1],list.get(i).replace(">", ""));
					
				}
			}
		}
		return countMap;
	}
	
	private static void addMiR(Map<String,StringBuilder> map, String key, String value){
		
		if(map.containsKey(key)){
			map.get(key).append(value+":");
		}
		else{
			StringBuilder sb = new StringBuilder();
			sb.append(value+":");
			map.put(key, sb);
		}
	}
	
	private static void add(Map<String,int[]> map, String key, int value){
		
		if(map.containsKey(key)){
			map.get(key)[0] += value;
		}
		else{
			int[] t = new int[1];
			t[0] = value;
			map.put(key, t);
		}
	}
	
	/**
	 * get only Drosha/Dicer duplexes back
	 */
	public static List<String> parseLax(List<String> list){

		List<String> back = new ArrayList<String>();
		for(int i = 0; i < list.size(); i++){
			if(list.get(i).startsWith(">")){
				if(list.size() > i+2){
					String[] f = list.get(i+2).split("\\s+")[0].split("&");
					int[] trailing1 = getTrailingNonBinding(f[0].toCharArray());
					int[] trailing2 = getTrailingNonBinding(f[1].toCharArray());
					if(trailing1[1] - trailing2[0] >= 1 && trailing1[1] - trailing2[0] <= 3 && trailing2[1] - trailing1[0] >= 1 && trailing2[1] - trailing1[0] <= 3){
						back.add(list.get(i));
						back.add(list.get(i+1));
						back.add(list.get(i+2));
					}
				}
			}
		}

		return back;
	}
	public static int[] getTrailingNonBinding(char[] sec){
		int[] back = new int[2];
		for(int i = 0; i < sec.length; i++){
			if(sec[i] == '.')
				back[0]++;
			else
				break;
		}
		for(int i = sec.length - 1; i > 0; i--){
			if(sec[i] == '.')
				back[1]++;
			else
				break;
		}	
		return back;
	}
	/**
	 * execute
	 */
	public static List<String> exec(String program, String parameter, String file){

		List<String> ll = new ArrayList<String>();
		String[] comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = program+" "+parameter+" < "+file;


		try {

			Runtime rt = Runtime.getRuntime();

			Process proc = null;

			proc = rt.exec(comand);
			InputStream stdin = proc.getInputStream();
			InputStreamReader isr = new InputStreamReader(stdin);
			BufferedReader br = new BufferedReader(isr);
			String line = null;
			while((line = br.readLine()) != null){
				ll.add(line);
			}
			proc.waitFor();
			br.close();
			isr.close();
			stdin.close();

		} catch (Throwable ee) {
			ee.printStackTrace();
		}


		return ll;
	}
	
	/**
	 *  Sort by Energy Ratio
	 * @param list
	 */
	
	public static void sort(List<RNAcofoldData> list){
		
		Comparator<RNAcofoldData> comparator = new Comparator<RNAcofoldData>(){
			public int compare(RNAcofoldData o1, RNAcofoldData o2) {
				if (o1.energyRatio > o2.energyRatio)
					return 1;
				if (o1.energyRatio < o2.energyRatio)
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);

	}
	/**
	 * Sort by energy
	 * @param list
	 */
	public static void sortEnergia(List<RNAcofoldData> list){
		
		Comparator<RNAcofoldData> comparator = new Comparator<RNAcofoldData>(){
			public int compare(RNAcofoldData o1, RNAcofoldData o2) {
				if (o1.energy > o2.energy)
					return 1;
				if (o1.energy < o2.energy)
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);

	}
	

	
	/**
	 * Aux. function for printing
	 * @param structureString
	 * @return
	 */
	public static int getPositionOfLastOpenBinding(String structureString){
		
		char[] asChar = structureString.toCharArray();
		for(int i = asChar.length - 1; i >= 0; i--){
			if(asChar[i] == '('){
				return i;
			}
		}
		return -1;
	}
	
	
	/**
	 * Aux. function for printing
	 * @param structureString
	 * @return
	 */
	public static int getPositionOfFirstCloseBinding(String structureString){
		
		char[] asChar = structureString.toCharArray();
		for(int i = 0; i < asChar.length ; i++){
			if(asChar[i] == ')'){
				return i;
			}
		}
		return -1;
	}
	/**
	 * 
	 * @param clust5p
	 * @param clust3p
	 * @param file
	 * @param extension5p
	 * @param maxNumber
	 */
	public  boolean representDuplexStacked(GFcluster clust5p, GFcluster clust3p, String file, int extension5p, int maxNumber){
		
		int spacer = 5;
		
		String[] ids = this.name.trim().split("#");
//		int id0 = Integer.parseInt(ids[0]);
		int id1 = Integer.parseInt(ids[1]);
//		int id2 = Integer.parseInt(ids[2]);
		int id3 = Integer.parseInt(ids[3]);
		
		int lineNumber = clust5p.seqs.size();
		if(clust3p.seqs.size() > lineNumber)
			lineNumber = clust3p.seqs.size();
		
		try {
			BufferedWriter writer = new BufferedWriter ( new FileWriter( file ));

			//			writer.write(Util.getVoidString(extension5p)+this.sec5p+Util.getVoidString(spacer)+
			//					Util.getVoidString(extension5p)+this.sec3p+"\n");


			for(int i = maxNumber; i >= 0; i--){
				if(clust5p.seqs.size() > i && i != id1){
					writer.write(Util.getVoidString(clust5p.start.get(i)+extension5p+spacer)+clust5p.seqs.get(i)+"\t"+clust5p.rc.get(i)+"\n");
				}
			}
			writer.write(Util.getVoidString(extension5p+spacer)+clust5p.prominentSeq+"\t"+clust5p.prominentRC+"\n");
			String structureString = null; // to align the last binding with the first binding of 3p
			
			// if mature is not dominant sequence
			if(id1 >= 0){
				writer.write(Util.getVoidString(clust5p.start.get(id1)+extension5p+spacer)+clust5p.seqs.get(id1)+"\t"+clust5p.rc.get(id1)+"\n\n");
				writer.write(Util.getVoidString(clust5p.start.get(id1)+extension5p+spacer)+this.struc5p+"\n");
				structureString = Util.getVoidString(clust5p.start.get(id1)+extension5p+spacer)+this.struc5p;
			}
			else{
				writer.write("\n"+Util.getVoidString(extension5p+spacer)+this.struc5p+"\n");
				structureString = Util.getVoidString(extension5p+spacer)+this.struc5p;
			}
			
			int positionOfLast = getPositionOfLastOpenBinding(structureString);
			
			if(positionOfLast < 0){
				IO.log(GVars.logFile, 4, "Error in structure "+this.name, true);
				writer.close();
				return false;
			}
			
			int open = getPositionOfFirstCloseBinding(this.struc3p);
			int voidLenRef = positionOfLast - (this.struc3p.length() - open) +1;
			if(voidLenRef < 0){
				IO.warning("Error in representing the structure "+this.name);
			}
			writer.write(Util.getVoidString(voidLenRef)+SeqUtil.getReverseSequence(this.struc3p)+"\n\n");

			int referenceStart = 0;
			if(id3 >= 0){
				referenceStart = clust3p.start.get(id3);
				writer.write(Util.getVoidString(voidLenRef)+SeqUtil.getReverseSequence(clust3p.seqs.get(id3))+"\t"+clust3p.rc.get(id3)+"\n");
				int voidLen = voidLenRef + clust3p.start.get(id3)+ (this.struc3p.length() - this.clust3p.prominentSeq.length());
				writer.write(Util.getVoidString(voidLen)+SeqUtil.getReverseSequence(clust3p.prominentSeq)+"\t"+clust3p.prominentRC+"\n");
			}
			else{
				writer.write(Util.getVoidString(voidLenRef)+SeqUtil.getReverseSequence(clust3p.prominentSeq)+"\t"+clust3p.prominentRC+"\n");
			}


			for(int i = 0; i <=maxNumber; i++){
				if(i == id3)
					continue;
				if(clust3p.seqs.size() > i ){
					// the start coordinates do not refer to the passenger strand but to the most frequent
					int voidLen = voidLenRef - (clust3p.start.get(i) -referenceStart )+ (this.struc3p.length() - this.clust3p.seqs.get(i).length());
					writer.write(Util.getVoidString(voidLen)+SeqUtil.getReverseSequence(clust3p.seqs.get(i))+"\t"+clust3p.rc.get(i)+"\n");
				}
			}


			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return true;
	}
	
	/**
	 * 
	 * @param list
	 */
	public static void sortBind(List<RNAcofoldData> list){
		
		Comparator<RNAcofoldData> comparator = new Comparator<RNAcofoldData>(){
			public int compare(RNAcofoldData o1, RNAcofoldData o2) {
				if (o1.bindings < o2.bindings)
					return 1;
				if (o1.bindings > o2.bindings)
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);
	}



}
