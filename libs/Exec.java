package libs;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


/**
 * This class should hold functions to launch external third party software <br>
 * Note: some programs might be launched from other classes, i.e. those that have there own class like RNAfold
 * @author Michael
 *
 */
public class Exec {

	
	public static boolean mv(String from, String to){
		
		
		String[] comand = null;
		comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = "mv "+from+" "+to;
		IO.log(GVars.logFile, 1, "mv "+from+" "+to, true);
		Runtime rt = Runtime.getRuntime();

		try{
			Process proc = rt.exec(comand);
			proc.waitFor();

		} catch (Throwable ee) {
			ee.printStackTrace();
			return false;
		}
		return true;
	}

	public static boolean cp(String from, String to){
		
		
		String[] comand = null;
		comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = "cp "+from+" "+to;
		IO.log(GVars.logFile, 1, "cp "+from+" "+to, true);
		Runtime rt = Runtime.getRuntime();

		try{
			Process proc = rt.exec(comand);
			proc.waitFor();

		} catch (Throwable ee) {
			ee.printStackTrace();
			return false;
		}
		return true;
	}
	
	
	public static boolean cmd(String command){
		
		
		String[] comand = null;
		comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = command;
		IO.log(GVars.logFile, 1,"will execute the following command: "+command, true);
		Runtime rt = Runtime.getRuntime();

		try{
			Process proc = rt.exec(comand);
			
			
			BufferedReader stdInput = new BufferedReader(new 
				     InputStreamReader(proc.getInputStream()));

				BufferedReader stdError = new BufferedReader(new 
				     InputStreamReader(proc.getErrorStream()));
//
				List<String> stdList = new ArrayList<String>(2);
//				List<String> errList = new ArrayList<String>(2);
				
				String s = null;
				// read the output from the command
//				System.out.println("Here is the standard output of the command:\n");
//
				while ((s = stdInput.readLine()) != null) {
				    System.out.println(s);
				}

				// read any errors from the attempted command
//				System.out.println("Here is the standard error of the command (if any):\n");
				while ((s = stdError.readLine()) != null) {
					stdList.add(s);
//				    System.out.println(s);
				}
				
				if(stdList.size() > 0 && !(stdList.get(0).startsWith("#"))){
					IO.warning("The command:"+comand[2]+" produced the following error output");
					for(String t : stdList){
						IO.writeToCommandLineL2(t);
					}
						
						
						
				}
				stdInput.close();
				stdError.close();
			proc.waitFor();

		} catch (Throwable ee) {
			ee.printStackTrace();
			return false;
		}
		return true;
	}

	/**
	 * 
	 * @param bowtieBinary
	 * @param file
	 * @param index
	 * @param overwrite
	 * @return: true --> index was built, false --> index was not built
	 */
	public static boolean makeBowtieIndex(String bowtieBinary, String file, String index, boolean overwrite){


		if(overwrite){

		}
		else{
			if(new File(index+".1.ebwt").exists() || new File(index+".1.ebwtl").exists() ){
				IO.log(GVars.logFile, 1, "Found index: "+index, true);
				return false;
			}
		}

		
		String[] comand = null;
		comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = bowtieBinary+" "+GVars.colorFlag+"  "+file+" "+index+GVars.colorIndex;
		IO.log(GVars.logFile, 1, "built command "+comand[2], true);
		Runtime rt = Runtime.getRuntime();

		try{
			Process proc = rt.exec(comand);
			proc.waitFor();

		} catch (Throwable ee) {
			ee.printStackTrace();
			return false;
		}
		return true;
	}


	
	
	
	public static String checkIfProgram(String execProc){
		
		if(System.getProperty("os.name").startsWith("Win")){
			Runtime rt = Runtime.getRuntime();
			String comand = execProc;
			try{


				Process proc = rt.exec(comand);
				String line = null;
				InputStream stdin = proc.getInputStream();
				InputStreamReader isr = new InputStreamReader(stdin);
				BufferedReader br = new BufferedReader(isr);
				
				while((line=br.readLine()) != null){

					System.out.println(line);
	
				}


				br.close();
				isr.close();
				stdin.close();
				proc.waitFor();
				return null;
			} catch (Throwable ee) {

			}
		}
		else{

			String[] comand = null;
			comand = new String[3];
			comand[0] = "/bin/bash";
			comand[1] = "-c";
			comand[2] = "which "+execProc;


			Runtime rt = Runtime.getRuntime();

			try{
				Process proc = rt.exec(comand);
				String line = null;
				InputStream stdin = proc.getInputStream();
				InputStreamReader isr = new InputStreamReader(stdin);

				BufferedReader br = new BufferedReader(isr);

				
				while((line=br.readLine()) != null){

					break;
					
				}
				
				br.close();
				isr.close();
				stdin.close();
				proc.waitFor();
				return line;

			} catch (Throwable ee) {
				ee.printStackTrace();
			}
		}
		return null;
	}
	
	/**
	 * launch the bowtie programm
	 * @param bowtieBinary
	 * @param indexPath
	 * @param inputFile
	 * @param parameters
	 * @param outFile
	 */
	public static void bowtieAlign(String bowtieBinary, String parameters){
		
	
		
		if(System.getProperty("os.name").startsWith("Win")){
			Runtime rt = Runtime.getRuntime();
			String comand = bowtieBinary+" "+parameters;
			try{


				Process proc = rt.exec(comand);
				String line = null;
				InputStream stdin = proc.getInputStream();
				InputStreamReader isr = new InputStreamReader(stdin);
				BufferedReader br = new BufferedReader(isr);

				
				while((line=br.readLine()) != null){

					System.out.println(line);
	
				}


				br.close();
				isr.close();
				stdin.close();
				proc.waitFor();
			} catch (Throwable ee) {

			}
		}
		else{

			String[] comand = null;
			comand = new String[3];
			comand[0] = "/bin/bash";
			comand[1] = "-c";
			comand[2] = bowtieBinary+"  "+parameters;


			Runtime rt = Runtime.getRuntime();

			try{
				Process proc = rt.exec(comand);
				

					BufferedReader stdError = new BufferedReader(new 
					     InputStreamReader(proc.getErrorStream()));
// 
					List<String> stdList = new ArrayList<String>(2);
					
					String s = null;
					while ((s = stdError.readLine()) != null) {
						stdList.add(s);
					}
					
					if(stdList.size() > 0 && !(stdList.get(0).startsWith("#"))){
						IO.warning("The command:"+comand[2]+" produced the following error output");
						for(String t : stdList){
							IO.writeToCommandLineL2(t);
						}
							
							
							
					}
					stdError.close();
				
				proc.waitFor();

			} catch (Throwable ee) {
				ee.printStackTrace();
			}
		}

	}
	
	
	
	
	/**
	 * launch the bowtie programm
	 * @param bowtieBinary
	 * @param indexPath
	 * @param inputFile
	 * @param parameters
	 * @param outFile
	 */
	public static void bowtieAlignDirect(String bowtieBinary, String parameters){
		
	
		
		if(System.getProperty("os.name").startsWith("Win")){
			Runtime rt = Runtime.getRuntime();
			String comand = bowtieBinary+" "+parameters;
			try{


				Process proc = rt.exec(comand);
				String line = null;
				InputStream stdin = proc.getInputStream();
				InputStreamReader isr = new InputStreamReader(stdin);
				BufferedReader br = new BufferedReader(isr);

				
				while((line=br.readLine()) != null){

					System.out.println(line);
	
				}


				br.close();
				isr.close();
				stdin.close();
				proc.waitFor();
			} catch (Throwable ee) {

			}
		}
		else{

			String[] comand = null;
			comand = new String[3];
			comand[0] = "/bin/bash";
			comand[1] = "-c";
			comand[2] = bowtieBinary+"  "+parameters;


			Runtime rt = Runtime.getRuntime();

			try{
				Process proc = rt.exec(comand);
				

				proc.waitFor();

			} catch (Throwable ee) {
				ee.printStackTrace();
			}
		}

	}

	/**
	 * 
	 * @param inFile
	 * @param parameter
	 * @param execFile
	 */
	public static void convertSRA(String inFile, String parameter ,String execFile){



		
		String[] comand = null;
		comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = execFile+"  "+parameter+" "+inFile;

		Runtime rt = Runtime.getRuntime();

		try{
			Process proc = rt.exec(comand);

			proc.waitFor();

		} catch (Throwable ee) {
			IO.warning("Error in SRA conversion (Exec.convertSRA () function). Will quit now");
			IO.log(GVars.logFile, 4, "Error in SRA conversion", true);
			System.exit(1);
//			ee.printStackTrace();
		}

	}
	
	
	


	/**
	 * 
	 * @param parameter
	 */
	public static void blast(String blastProgram, String parameter  ){

		Runtime rt = null;
		Process proc = null;
		String comand1 = blastProgram+" "+parameter ;
		String[] comand = new String[3];
		comand[0] = "/bin/bash";
		comand[1] = "-c";
		comand[2] = comand1;
		System.out.println(comand1);
		rt = Runtime.getRuntime();
		
		try {
//			proc = pb.start();
			proc = rt.exec(comand);
			proc.waitFor();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public static void removeDir(String dir){

		Runtime rt = null;
		Process proc = null;
		String comand = "rm -r "+dir+"";

		System.out.println(comand);
		

		rt = Runtime.getRuntime();
		String[] command = new String[3];
		command[0] = "/bin/bash";
		command[1] = "-c";
		command[2] = comand;
		
		try {

			proc = rt.exec(command);
			proc.waitFor();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static String wget(String url, String parameters, String pwd){
		
		Runtime rt = null;
		Process proc = null;
		
		String outfile = pwd + File.separator + new File(url).getName();
		if(new File(outfile).exists()){
			IO.writeToCommandLineL1(outfile+" exists");
			return outfile;
		}
		String comand = "cd "+pwd+"; wget "+parameters+" "+url+" ";

		String[] command = new String[3];
		command[0] = "/bin/bash";
		command[1] = "-c";
		command[2] = comand;
		
		System.out.println(comand);
		
		rt = Runtime.getRuntime();
		
		try {
//			proc = pb.start();
			proc = rt.exec(command);
			int i = proc.waitFor();
			
			return pwd + File.separator + new File(url).getName();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	

}
