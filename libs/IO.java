package libs;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;



public class IO {

	// standard output parameters
	static public int outWidth = 70; // the width of the output


	////////////////////////////////////////////////////////////////////////////////////////////////
	//////////  	OUTPUT TO THE STDOUT
	///////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * write out warning to the standard output
	 * @param line
	 */
	public static void warning(String line){
		System.out.println("\n"+Util.getCharString(outWidth, '#'));
		System.out.println("\n     "+line+"     \n");
		System.out.println("\n"+Util.getCharString(outWidth, '#'));
	}


	public static void writeToCommandLineL1(String line){

		System.out.println("\n             "+line+"\n");

	}
	public static void writeToCommandLineL2(String line){

		System.out.println("               "+line+"");

	}

	public static void writeToCommandLineBlockOpen(String line){

		System.out.println(Util.getCharString(outWidth, '-'));
		System.out.println(Util.getCharString(10, '-')+Util.getVoidString(outWidth-20)+Util.getCharString(10, '-'));
		System.out.println("           "+line+Util.getVoidString(60-(line.length() + 5)));

	}

	public static void writeToCommandLineBlockClose(String line){

		System.out.println("\n           "+line);
		System.out.println(Util.getCharString(10, '-')+Util.getVoidString(outWidth-20)+Util.getCharString(10, '-'));
		System.out.println(Util.getCharString(outWidth, '-')+"\n");

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/**
	 *  Reads the config file and writes the parameters into the info hash
	 * @param file
	 */
	public static Map<String,List<String>> readConfigFile(String file){

		Map<String,List<String>> back = new Hashtable<String,List<String>>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(file));
			String line = null;
			while ((line = in.readLine()) != null) {
				if(line.contains("=")){
					String[] f = line.trim().split("=",2);
					if(back.containsKey(f[0])){
						back.get(f[0]).add(f[1]);
					}
					else{
						List<String> temp = new ArrayList<String>(2);
						if(f.length > 1)
							temp.add(f[1]);
						else
							temp.add(" ");
						back.put(f[0], temp);
					}
				}
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.warning("Config file not found! Will quit now.");
			System.exit(1);
		}
		return back;

	}


	/**
	 * Using stream to copy file
	 * @param source
	 * @param dest
	 * @throws IOException
	 */
    public static void copyFileUsingStream(File source, File dest) throws IOException {
        InputStream is = null;
        OutputStream os = null;
        try {
            is = new FileInputStream(source);
            os = new FileOutputStream(dest);
            byte[] buffer = new byte[1024];
            int length;
            while ((length = is.read(buffer)) > 0) {
                os.write(buffer, 0, length);
            }
        } finally {
            is.close();
            os.close();
        }
    }
    
    
    /**
     * copy a gz file and gunzip it
     * @param input
     * @param output
     */
    public static boolean copyUnzipFile(String input, String output){

    	
    	if(new File(output).exists()){
    		IO.writeToCommandLineL1(output+" exists");
    		return true;
    	}
    	
    	try{
    		BufferedWriter writer = new BufferedWriter (new FileWriter (output));
    		BufferedReader reader;
    		if(input.endsWith("gz")){
    			reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(input) ) ) );
    		}
    		else{
    			reader = new BufferedReader (new FileReader(input));
    		}
    		String line = null;

    		while((line = reader.readLine()) != null){
    			writer.write(line+"\n");
    		}
    		writer.close();
    		reader.close();
    		return true;
    	}catch(FileNotFoundException e){
    		warning(input+" not found (IO.manipulateFastaFile() function). Will quit now!");
    		log(GVars.logFile, 4, "File not found: "+new File(input).getName(), true);
    		System.exit(1);
    	} catch (IOException e) {
    		// TODO Auto-generated catch block
    		warning("Problems with file: "+input+".");
    		log(GVars.logFile, 3, "Probles with file: "+new File(input).getName()+".  ", true);
    		e.printStackTrace();
    	}
    	return false;
    }

	/**
	 * Copy a zip file
	 * @param zipFile
	 * @param newFile
	 * @throws IOException
	 */
	public static void copyZipFile(File zipFile, File newFile) {
		


		try {
			ZipFile zipSrc = new ZipFile(zipFile);
			ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(newFile));
			Enumeration srcEntries = zipSrc.entries();
			while (srcEntries.hasMoreElements()) {
				ZipEntry entry = (ZipEntry) srcEntries.nextElement();
				ZipEntry newEntry = new ZipEntry(entry.getName());

				zos.putNextEntry(newEntry);
				BufferedInputStream bis = new BufferedInputStream(zipSrc
						.getInputStream(entry));

				while (bis.available() > 0) {
					zos.write(bis.read());
				}
				zos.closeEntry();

				bis.close();


			}
			zos.finish();
			zos.close();
			zipSrc.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile, 4, "Error copying a zip file in IO.copyZipFile. Will quit now", true);
			System.exit(1);
		}

	}
	
	/**
	 * copies a file
	 * @param in
	 * @param out
	 * @param append --> true: append to an existing file
	 */
	public static void copy(String in, String out, boolean append){

		if(new File(in).exists()){
			try {
				BufferedWriter writer = new BufferedWriter (new FileWriter(out, append));
				BufferedReader reader = new BufferedReader(new FileReader(in));

				String line = null;
				while((line = reader.readLine()) != null){
					writer.write(line+"\n");
				}
				reader.close();
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * copies a file
	 * @param in
	 * @param out
	 * @param append --> true: append to an existing file
	 */
	public static void copy(String in, String out, boolean append, boolean copyHeader){

		if(new File(in).exists()){
			try {
				BufferedWriter writer = new BufferedWriter (new FileWriter(out, append));
				BufferedReader reader = new BufferedReader(new FileReader(in));

				String line = null;
				if(!(copyHeader)){
					reader.readLine();
				}
				while((line = reader.readLine()) != null){
					writer.write(line+"\n");
				}
				reader.close();
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public static void manipulateFastaFile(String input, String output, String add ){
		
		try {

			BufferedWriter writer  = new BufferedWriter(new FileWriter(output));
			BufferedReader reader;
			if(input.endsWith("gz")){
				reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(input) ) ) );
			}
			else{
				reader = new BufferedReader (new FileReader(input));
			}
			String line = null;

			while((line = reader.readLine()) != null){

				if(line.startsWith(">")){
					String[] f = line.split("\\s+");
					String o = f[0]+add;
					writer.write(o+"\n");
				}
				else{
					writer.write(line+"\n");
				}

			}
			writer.close();
			reader.close();
		}catch(FileNotFoundException e){
			warning(input+" not found (IO.manipulateFastaFile() function). Will quit now!");
			log(GVars.logFile, 4, "File not found: "+new File(input).getName(), true);
			System.exit(1);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			warning("Problems with file: "+input+".");
			log(GVars.logFile, 3, "Probles with file: "+new File(input).getName()+".  ", true);
			e.printStackTrace();
		}
	}
	
	/**
	 * Make a subfile out of the fastq input file 
	 * @param file: input file. must be fastq
	 * @param outFile: output file. will be fastq
	 * @param entries: number of reads in output file
	 */
	public static void makeSubFastq(String file, String outFile, int entries){
		
		
		String[] files = file.split(":");
		
		BufferedReader reader;

		try {

			BufferedWriter writer  = new BufferedWriter(new FileWriter(outFile));
			int c = 0;
			for(String file1 : files){
//				reader = new BufferedReader(new FileReader(file1));
				System.out.println();
				writeToCommandLineL2("Reading file: "+file1);
				if(file1.endsWith("gz")){
					 reader = new BufferedReader(new InputStreamReader(	new GZIPInputStream( new FileInputStream(file1) ) ) );
				}
				else{
					 reader = new BufferedReader (new FileReader(file1));
				}
				String line = null;

				while((line = reader.readLine()) != null){

					String seq = reader.readLine();
					if(seq == null){
						warning("Read Sequence expected after line: "+line+" but not found!!!");
						break;
					}
					seq = seq.trim();

//					System.out.println(seq);
					String idDumm = reader.readLine();
					String qual = reader.readLine();
					if(c < entries){
						writer.write(line+"\n");
						writer.write(seq+"\n");
						writer.write(idDumm+"\n");
						writer.write(qual+"\n");
					}
					else{
						break;
					}
					c++;
				}
				reader.close();
			}
			writer.close();
		}catch(FileNotFoundException e){
			warning(file+" not found (IO.makeSubFastq() function). Will quit now!");
			log(GVars.logFile, 4, "File not found: "+new File(file).getName(), true);
			System.exit(1);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			warning("Problems with file: "+file+". Very likely the number of lines is not a multiple of 4 (IO.makeSubFastq() function).");
			log(GVars.logFile, 3, "Probles with file: "+new File(file).getName()+". Very likely the number of lines is not a multiple of 4", true);
			e.printStackTrace();
		}
	}






	/**
	 * write into the log file
	 * @param file
	 * @param label: 1 --> INFO; 2 --> SUCCESS; 3 --> WARNING, 4 --> ERROR; 5 --> Backvalues
	 */
	public static void log(String file,int label, String message, boolean append){

		if(file != null){

			if(label == 1){
				Write.writeString(file, new Timestamp(new java.util.Date().getTime())+" INFO: "+  message, append);
			}
			else if (label == 2){
				Write.writeString(file, new Timestamp(new java.util.Date().getTime())+" SUCCESS: "+  message, append);
			}
			else if (label == 3){
				Write.writeString(file, new Timestamp(new java.util.Date().getTime())+" WARNING: "+  message, append);
			}
			else if (label == 4){
				Write.writeString(file, new Timestamp(new java.util.Date().getTime())+" ERROR: "+  message, append);
			}
			else if (label == 5){
				Write.writeString(file, new Timestamp(new java.util.Date().getTime())+" BACKVALUE: "+  message, append);
			}
			else if (label == 6){
				Write.writeString(file, new Timestamp(new java.util.Date().getTime())+" OTHER: "+  message, append);
			}

		}
	}
	
	static public boolean saveUrl(final String filename, final String urlString){
	   
		if(new File(filename).exists()){
			IO.writeToCommandLineL1(filename+" exists already");
			return true;
		}
		BufferedInputStream in = null;
		FileOutputStream fout = null;
		try {
			try {
				in = new BufferedInputStream(new URL(urlString).openStream());
				fout = new FileOutputStream(filename);

				final byte data[] = new byte[1024];
				int count;
				while ((count = in.read(data, 0, 1024)) != -1) {
					fout.write(data, 0, count);
				}
			} finally {
				if (in != null) {
					in.close();
				}
				if (fout != null) {
					fout.close();
				}
			}
			return true;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IO.log(GVars.logFile,4, "problem with URL: "+urlString, true);
			IO.warning("problem with URL: "+urlString);
//			e.printStackTrace();
		}
		return false;
	}
	
	
	/**
	 * the funtion splits each line by '\t' and writes out into a separate file all lines which start at column 'col' with 'searchPat'
	 * @param infile
	 * @param outfile
	 * @param searchPath
	 * @param col
	 */
	public static void parseOutSubfile(String infile, String outfile, String searchPat, int col, boolean header){

		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter(outfile));
			BufferedReader reader = new BufferedReader (new FileReader(infile));
			String line = null;
			String headStr=null;

			if(header){
				headStr = reader.readLine();
				writer.write(headStr+"\n");
			}
			
			while((line = reader.readLine()) != null){
				String[] f = line.split("\t");
				if(f[col].startsWith(searchPat)){
					writer.write(line+"\n");
				}
			}
			
			reader.close();
			writer.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
//	public static void getFTP(){
//		
//		
//		private static final int BUFFER_SIZE = 4096 ;			
//		long startTime = System.currentTimeMillis() ; 
//		
//		String ftpUrl = ftp://**username**:**password**@filePath ; 
//		String file= "filename" ; // name of the file which has to be download 
//		String host = host_name ; //ftp server 
//		String user = "username" ; //user name of the ftp server 
//		String pass = "**password" ; // password of the ftp server 
//		String savePath = "c:\\" ; 
//		ftpUrl = String.format(ftpUrl, user, pass, host) ; 
//		System.out.println("Connecting to FTP server") ; 
//		
//		try{
//			URL url = new URL(ftpUrl) ;
//			URLConnection conn = url.openConnection() ; 
//			InputStream inputStream = conn.getInputStream() ; 
//			long filesize = conn.getContentLength() ; 
//			System.out.println("Size of the file to download in kb is:-" + filesize/1024 ) ;
//			FileOutputStream outputStream = new FileOutputStream(savePath) ; 
//			byte[] buffer = new byte[BUFFER_SIZE] ; 
//			int bytesRead = -1 ; 
//			while ((bytesRead = inputStream.read(buffer)) != -1) 
//			{
//				outputStream.write(buffer, 0, bytesRead) ; } 
//			long endTime = System.currentTimeMillis() ; 
//			System.out.println("File downloaded") ; 
//			System.out.println("Download time in sec. is:-" + (endTime-startTime)/1000) ; 
//			outputStream.close() ; 
//			inputStream.close() ; }
//		catch (IOException ex){ ex.printStackTrace() ; } } 
//
//			Read more: http://mrbool.com/java-ftp-how-to-download-file-with-java/29831#ixzz3aQhArsya
//	}
//}


	public static void noAccess(String logFile, String outFile){
		IO.log(logFile, 4, "Could not access "+new File(outFile).getName(), true);
		IO.warning( "Could not access "+outFile);
	}
	public static void noFile(String logFile, String outFile){
		IO.log(logFile, 4, "File not found: "+new File(outFile).getName(), true);
		IO.warning( "File not found: "+outFile);
	}

	public static void jobInterupt(String logFile, String launchCommand){
		IO.log(logFile, 4, "the following launch command failed: "+launchCommand, true);
		IO.warning( "the following launch command failed: "+launchCommand);
	}
}
