package miRNAgFree;

import java.io.File;
import java.util.List;

import libs.GVars;
import libs.IO;
import libs.Preproc;
import libs.Stat;



public class MiRNAgFree {

	public static void main(String[] args) {
		// TODO Auto-generated method stub

		Vars.welcome();
		////////////////////////////////////////////////////////////////////////
		//////////////
		/////////////      1) reads the parameters either from the command line or from a config file
		/////////////      2) sets the default values for those that do not appear
		/////////////      3) prepare the output directories	
		Vars.getParameters(args);
		
		///////////////////////////////////////////////////////////////////////////////
		///////     The input is not previously aligned (SAM/BAM) (i.e. if doAlignment==true)
		//////        --> several data preparation steps & preprocessing of the reads
		//////      * convert SRA to Fastq
		//////      * prepare SOLID data
		/////       * convert the input file to an internal sRNAbench format fasta file
		/////       * trimm the adapter (if adapter= is given), removes barcodes and 3' end random oligos (ligase bias reduction protocols)
		/////	
		boolean preprocessed = Preproc.preprocessing();
//		Results.readLengthFile = GVars.stat+File.separator+"readLengthAnalysis.txt";
		///////////////////////////////////////////////////////////////////////////////
		///////     Make length distribution of adapter trimmed reads
		//////      
		if(preprocessed)
			Stat.makeReadLengthStat(GVars.origInput, GVars.input);
		
		////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////
		//////////////
		/////    prepare the analysis
		////     i) read microRNA data for those that are used as putative homologous
		////     ii) set the thresholds as a function of the 'mode' (i.e. strict or lax)
		GfreeHelper.prepare();

		//////////////////
		//// apply the filter 
		Preproc.applyFilter(Vars.libsFilter, GVars.input, GVars.output, false);

		///////////////////
		/// cluster the reads --> each GFcluster object holds a read cluster
		List<GFcluster> clusterList = Cluster.cluster(GVars.output+File.separator+"reads.fa");
		
		IO.writeToCommandLineL1("Found "+ clusterList.size()+" clusters");
		
		///////////////////////////
		///  1) apply thresholds to each putative guide cluster
		///  2) determine the most likely passenger strand
		List<GFresult> resList = DetectMiRNA.detectMicroRNAs(clusterList);
		GFresult.sortBiggerToSmaller(resList);
		IO.writeToCommandLineL1("Found "+ resList.size()+" miRNAs");
		
		//////////////////////////
		/// write out the fasta files, expression files and graphical representation
		GfreeHelper.writeOutResults(resList, Preproc.reads);
		
	}

}
