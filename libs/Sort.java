package libs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class Sort {

	
	public double sortValue;
	public String line;
	public Sort(double sortValue, String line){
		
		this.sortValue = sortValue;
		this.line = line;
	}

	/**
	 * 
	 * @param file
	 * @param column --> 0 based!!
	 * @param header
	 * @return
	 */
	public static List<Sort> sortListBigger2Smaller(String file, int column, boolean header){
		
		List<Sort> back = new ArrayList<Sort>();
		String head = null;
		try{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			if(header){
				head = reader.readLine();
			}
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				Sort s = new Sort(Double.parseDouble(f[column]),line);
				back.add(s);
			}
			reader.close();
			sortBiggerToSmaller(back);
			writeOut(file, back, head);
			
		}
		catch (IOException e){
//			e.printStackTrace();
		}
		return back;
	}

	
	
	/**
	 * 
	 * @param file
	 * @param column --> 0 based!!
	 * @param header
	 * @return
	 */
	public static List<Sort> sortListSmaller2Bigger(String file, int column, boolean header){
		
		List<Sort> back = new ArrayList<Sort>();
		String head = null;
		try{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			if(header){
				head = reader.readLine();
			}
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				Sort s = new Sort(Double.parseDouble(f[column]),line);
				back.add(s);
			}
			reader.close();
			sortSmallerToBigger(back);
			writeOut(file, back, head);
		}
		catch (IOException e){
			e.printStackTrace();
		}
		return back;
	}

	
	
	public static void sortSmallerToBigger(List<Sort> list){


		Comparator<Sort> comparator = new Comparator<Sort>(){
			public int compare(Sort o1, Sort o2) {
				if (o1.sortValue > o2.sortValue)
					return 1;
				if (o1.sortValue < o2.sortValue)
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);

	}

	public static void sortBiggerToSmaller(List<Sort> list){
		
		
		Comparator<Sort> comparator = new Comparator<Sort>(){
			public int compare(Sort o1, Sort o2) {
				if (o1.sortValue < o2.sortValue)
					return 1;
				if (o1.sortValue > o2.sortValue)
					return -1;
				return 0;
			}
		};
		Collections.sort(list, comparator);

	}

	
	public static void writeOut(String file, List<Sort> list, String header){
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			if(header != null)
				writer.write(header+"\n");
			for(Sort sort : list){
				writer.write(sort.line+"\n");
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	/**
	 * 
	 * @param list
	 * @param column: --> 0-based!!!
	 * @return
	 */
	public static List<Sort> getSortList(List<String> list, int column){
		
		List<Sort> back = new ArrayList<Sort>();
		for(String str : list){
			String[] f = str.split("\t");
			back.add(new Sort(Double.parseDouble(f[column]), str));			
		}
		return back;
	}
	
	/**
	 * gives back the original list sorted
	 * @param list
	 * @return
	 */
	public static List<String> getOriginalList(List<Sort> list){
		
		List<String> back = new ArrayList<String>();
		for(Sort str : list){
			back.add(str.line);
		}
		return back;
	}
}
