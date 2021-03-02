package watseqNAM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;



public class ExtractNucmerSNPs_parallel {

	HashMap<String, Vector<int[]>> nucmerAlignments;  //<chr, <{alignmentStart,alignmentEnd}>
	HashMap<String, HashMap<Integer, String[]>> nucmerSnps; //<chr, <pos, {csAllele, paragonAllele}>>
	
	HashMap<String, HashMap<String, String>> axiomData;   // <markerName, <recombinant, haplotype>>
	HashMap<String, Vector<String>> populations;  //<watseqAccession, <recombinant>>
	
	HashMap<String, HashMap<Integer, String>> axiomMarkerPositions;  //<chr, <pos,markername>>
	HashMap<String, HashMap<String, Vector<Integer>>> availablePositions; //<chr,<accession, <pos>>>
	
	
	HashMap<String, BufferedWriter> outputVCFs;
	
	
	
	
	public static void main(String[] args) {
		

		
		
		CLI cli = new CLI();
		cli.parseOptions(args);
		
		if( !cli.hasOption("i")|| !cli.hasOption("o") || !cli.hasOption("s") || !cli.hasOption("a") || !cli.hasOption("m")  || !cli.hasOption("c")) {
			System.out.println("java -jar WatSeq_VCFsForNAM.jar -i input.vcf_files -o outputDirectory -s nucmerSNPs.txt -a nucmerAlignments.txt -m axiomMap.csv_files -c nameConversionTable.txt -p outputFilePrefix");
			
		}else {
			File outputFile = new File(cli.getArg("o"));
			
			Vector<File> vcfs = new Vector<File>();
			
			for( Iterator<String> iterator = cli.getArgs("i").iterator(); iterator.hasNext();) {
				vcfs.add(new File(iterator.next()));
			}
			
			
			Vector<File> nucmerSNPFiles = new Vector<File>();
			for(Iterator<String> iterator = cli.getArgs("s").iterator(); iterator.hasNext();) {
				nucmerSNPFiles.add(new File(iterator.next()));
			}
			
			Vector<File>nucmerAlignmentFiles = new Vector<File>();
			for(Iterator<String> iterator = cli.getArgs("a").iterator(); iterator.hasNext();) {
				nucmerAlignmentFiles.add(new File(iterator.next()));
			}
			
			Vector<File> axiomFiles = new Vector<File>();
			for(Iterator<String> iterator = cli.getArgs("m").iterator(); iterator.hasNext();) {
				axiomFiles.add(new File(iterator.next()));
			}
			
			
			File nameConversionTableFile = new File(cli.getArg("c"));
			String outputFilePrefix = "watseq_test";
			
			if(cli.hasOption("p")) {
				outputFilePrefix = cli.getArg("p");
			}
			
			try {
				//mergeFiles(nucmerSNPFile, nucmerAlignmentFile, vcfInput, vcfOutput);
				ExtractNucmerSNPs_parallel e = new ExtractNucmerSNPs_parallel(nucmerAlignmentFiles, nucmerSNPFiles, axiomFiles,  nameConversionTableFile);
				e.imputeSNPs(vcfs, outputFile, outputFilePrefix);
			} catch (Exception e) {
				e.printStackTrace();
			}

		}
			
		
	}

	public ExtractNucmerSNPs_parallel(Vector<File> nucmerAlignmentFiles, Vector<File> nucmerSNPFiles, Vector<File> axiomFiles, File nameConversionTableFile) throws IOException{
		
		
		
		readAxiomFiles(axiomFiles,nameConversionTableFile);
		
		
		
		this.outputVCFs = new HashMap<String, BufferedWriter>();
		
		
		
		for(Iterator<File> iterator = nucmerAlignmentFiles.iterator(); iterator.hasNext();) {
			File nucmerAlignmentFile = iterator.next();
			readNucmerAlignments(nucmerAlignmentFile);
			System.out.println("read " +  nucmerAlignmentFile.getName());
		}
		
		for(Iterator<File> iterator = nucmerSNPFiles.iterator(); iterator.hasNext();) {
			File nucmerSNPFile = iterator.next();
			readNucmerSNPs(nucmerSNPFile);
			System.out.println("read " +  nucmerSNPFile.getName());
		}
			
		
		
		
		
	}
	
	
	
	
	public void imputeSNPs(Vector<File> inputVCFs, File outputDir, String outputFilePrefix)throws IOException{
		
		if(!outputDir.exists()) {
			outputDir.mkdirs();
		}
		
		
		for(Iterator<String> iterator = populations.keySet().iterator(); iterator.hasNext();) {
			String accession = iterator.next();
			
			BufferedWriter outVCF = new BufferedWriter(new FileWriter(new File( outputDir, accession + "_" + outputFilePrefix + ".vcf" )));
			this.outputVCFs.put(accession, outVCF);
			this.outputVCFs.get(accession).write("##fileformat=VCFv4.2");
			this.outputVCFs.get(accession).newLine();
			this.outputVCFs.get(accession).write("##INFO=<ID=axiomMarker,Number=1,Type=String,Description=\"imputation marker, reference Paragon\">");
			this.outputVCFs.get(accession).newLine();
			this.outputVCFs.get(accession).write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
			
			
			this.outputVCFs.get(accession).newLine();
			this.outputVCFs.get(accession).write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			
		
			
			for(Iterator<String> iterator2 = this.populations.get(accession).iterator(); iterator2.hasNext();) {
				String recombinant = iterator2.next();
				
				this.outputVCFs.get(accession).write("\t" + recombinant);
			}
			this.outputVCFs.get(accession).newLine();
			
		}
		
		
		
		for(Iterator<File> vcfIterator = inputVCFs.iterator(); vcfIterator.hasNext(); ) { //iterate through vcfs, i.e. read one chromosome at a time.
			
			
			
			File inputVCF = vcfIterator.next();
			
			BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream (inputVCF))));
			
			HashMap<String, Integer> columnNumbers = new HashMap<String, Integer>();
			 
			
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if( inputline.startsWith("##")) {
					continue;
				}
				String[] split = inputline.split("\t");
				
				
				
				//record the column headers. this is to know from which column we need to extract the genotype
				
				if(inputline.startsWith("#CHROM")) {
					for( int i = 0; i< split.length; i++) {
						if( populations.containsKey(split[i])) {
							columnNumbers.put(split[i], i) ;
						}
					}
					continue;
				}
				
				
				
				String chr = split[0].substring(3);
				int pos = Integer.parseInt(split[1]);
				String snp_id = split[2];
				
				//see if we have a paragon allele.
				String paragonAllele = this.getParagonAllele(chr, pos, split[3]);
								
				if( paragonAllele == null) {
					continue;
				}
				
				//get the genotype for our accession.
				String gt;   
				int genotypeFieldIndex = 0;
				String[] formatSplit = split[8].split(":");
				for( int i = 0; i< formatSplit.length; i++) {
					if(formatSplit[i].equalsIgnoreCase("GT")) {
						genotypeFieldIndex = i;
						break;
					}
				}
				
				
				String vcfRefAllele = split[3];
				String vcfAltAllele = split[4];
				
				
				// skip SNPs where the paragon allele is either ref nor alt from the vcf -> keep it bi-allelic
				if( !paragonAllele.equalsIgnoreCase(vcfRefAllele) && ! paragonAllele.equalsIgnoreCase(vcfAltAllele)) {
					continue;
				}
				
				
				//in the output the paragon allele should be the ref allele. So the gt needs to be switched if paragon was alt. 
				boolean paragonWasAlt = false;
				if(paragonAllele.equalsIgnoreCase(vcfAltAllele)) {
					paragonWasAlt = true;
				}
				
			
				
				
				
				// for each accession get watseq allele
					
				for( Iterator<String> iterator1 = this.populations.keySet().iterator(); iterator1.hasNext();) {
					String accession = iterator1.next();
					
					
					if(!columnNumbers.containsKey(accession)) {
						continue;
					}
					
					gt = split[columnNumbers.get(accession)].split(":")[genotypeFieldIndex];
					boolean watSeqIsHet = false;
					
					if(gt.equalsIgnoreCase("1/0") || gt.equalsIgnoreCase("1|0")||gt.equalsIgnoreCase("0/1") || gt.equalsIgnoreCase("0|1") ) {
						watSeqIsHet = true;
					}
					
					if( paragonWasAlt ) {
						
						if( gt.equalsIgnoreCase("1/1") || gt.equalsIgnoreCase("1|1") ) {
							gt = "0/0";
						}
						if( gt.equalsIgnoreCase("0/0") || gt.equalsIgnoreCase("0|0") ) {
							gt = "1/1";
						}
						
					}
					
					
					//find the next axiomMarker to this position
					String imputationMarker = this.getImputationMarker(chr, pos, accession);
					
					
					
					
					//this.outputVCFs.get(accession).write(markerName + "\t" + alleles+"\t" + chr + "\t" + pos +"\t+\tRefSeq-v1.0\tJIC\tNA\tNA\tJICcollection\tNA");
					this.outputVCFs.get(accession).write(chr + "\t" + pos + "\t" + snp_id + "\t");
					this.outputVCFs.get(accession).write(paragonAllele + "\t");
					
					String outputAltAllele = vcfAltAllele;
					if(paragonWasAlt) {
						outputAltAllele = vcfRefAllele;
					}
					this.outputVCFs.get(accession).write(outputAltAllele + "\t");
					this.outputVCFs.get(accession).write(".\t.\taxiomMarker=" + imputationMarker + "\t");
					this.outputVCFs.get(accession).write("GT");		
					
					
					
					for(Iterator<String> iterator2 = this.populations.get(accession).iterator(); iterator2.hasNext();) {
						String recombinant = iterator2.next();
						String accessionAllele = null;
						
						try {
							accessionAllele = this.axiomData.get(imputationMarker).get(recombinant);  //there will be a nullpointer if there is no marker on this chromosome for this accession.
						}catch (NullPointerException e) {}
						
						
						String outputGT = "0/0";
						
						
						if(accessionAllele==null) {
							outputGT = ".";
						}else if (accessionAllele.equalsIgnoreCase("H") || accessionAllele.equalsIgnoreCase("X") || watSeqIsHet) {
							outputGT="0/1";
						}else if (accessionAllele.equalsIgnoreCase("B")) {
							outputGT = "1/1";
						}
						
						this.outputVCFs.get(accession).write("\t" + outputGT);
					}
					
					this.outputVCFs.get(accession).newLine();
					
				}
				
				
				
				
				
				
				
				
				
			}

			in.close();
			
			System.out.println("read " + inputVCF.getName());	
		}
		
		
		for(Iterator<String> iterator = this.outputVCFs.keySet().iterator(); iterator.hasNext();) {
			this.outputVCFs.get(iterator.next()).close();
		}
		
		
		
		
	}
	
	
	/**
	 * 
	 * extract the paragon allele from the recorded nucmer data.
	 * 
	 * If a SNP was recorded the alternative paragon allele is returned.
	 * If there was no SNP but the position is in an aligned region, the Chinese Spring allele is returned.
	 * If the position is not in an aligned region it returnes null.
	 * It also returns null for a SNP where the refAllele from the input to this method is not consistent to the ref allele from nucmer.
	 * 
	 * @param chr chromosome of the snp
	 * @param pos position of the snp
	 * @param refAllele the reference allele.
	 * @return the paragon allele.
	 */
	private String getParagonAllele(String chr, int pos, String refAllele) {
		
		
		if(this.nucmerSnps.containsKey(chr) && this.nucmerSnps.get(chr).containsKey(pos)) {
			String[] allele = this.nucmerSnps.get(chr).get(pos);
			
			if(allele[0].equalsIgnoreCase(refAllele)) {
				return allele[1];
			}else {
				System.out.println("[debug] Not matching reference allele for " + chr + ":" + pos + "\t" + refAllele + "\t" + allele[0] );
				return null;
			}
			
			
		}else {
			
			for(Iterator<int[]> iterator = this.nucmerAlignments.get(chr).iterator(); iterator.hasNext();) {
				int[] a = iterator.next();
				if( pos >= a[0] && pos <= a[1]) {
					return refAllele;
				}else {
					if( a[0] > pos) {
						break;
					}
				}
				
			}
			
			
		}
		
		return null;
	}
	
	private String getImputationMarker( String chr, int pos , String accession) {
		String imputation = "-";
		
		Vector<Integer> v = this.availablePositions.get(chr).get(accession);
		if (v == null) {
			return imputation;
		}
		int[] a = {0,v.size()-1} ;
		while(a[1]-a[0]>1) {
			int index = a[0] + (a[1] -a[0])/2;
			if( v.get(index) > pos ) {
				a[1] = index;
			}else {
				a[0] = index;
			}
		}
		
		//System.out.println(a[0] + "\t" + a[1] + "\t"+ v.get(a[0]) + "\t" + pos + "\t" + v.get(a[1]));
		
		int markerPosition = v.get(a[0]);
		if( v.get(a[1]) - pos < pos - v.get(a[0])) {
			markerPosition = v.get(a[1]);
		}
		
		imputation = axiomMarkerPositions.get(chr).get(markerPosition);
		//System.out.println(imputation);
		return imputation;
	}
	
	
	/**
	 * 
	 * read the data from the axiom files
	 * 
	 * @param axiomFiles
	 * 		a list of the axiom files
	 * @param nameConversionTable
	 * 		the table with the name conversions that translates the prefix of the axiom file to the accession name that can be searched for in the vcf.
	 * 		The prefix is in column 0 and the accession is in column 4
	 * @throws IOException
	 */
	private void readAxiomFiles(Vector<File> axiomFiles, File nameConversionTableFile)throws IOException{
		
		HashMap<String, String> nameConversionTable = this.readNameConversionTable(nameConversionTableFile);
		
		if(this.outputVCFs == null) {
			this.outputVCFs = new HashMap<String, BufferedWriter>();
		}
		
		if(this.axiomData == null) {
			this.axiomData = new HashMap<String, HashMap<String, String>>();
		}
		
		if( this.axiomMarkerPositions == null) {
			this.axiomMarkerPositions = new HashMap<String, HashMap<Integer, String>>();
		}
		
		if(this.availablePositions == null) {
			this.availablePositions = new HashMap<String, HashMap<String,Vector<Integer>>>();
		}
		
		if( this.populations == null) {
			this.populations = new HashMap<String, Vector<String>>();
		}
		
		for(Iterator<File> iterator = axiomFiles.iterator(); iterator.hasNext();) {
			File axiomFile = iterator.next();
			String accession = nameConversionTable.get(axiomFile.getName().split("_")[0]); //this is the watseq parent
			
			populations.put(accession, new Vector<String>());
			
				
			BufferedReader in = new BufferedReader(new FileReader(axiomFile));
			
			String[] header = in.readLine().split(",");
			
			
			for( int i = 5; i < header.length; i++) {
				String recombinant = header[i];
			
				this.populations.get(accession).add(recombinant); //add the order of the recombinants
				//System.out.println("Populations: adding " + recombinant + " to " + accession);
			}
			
			
			
			
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				
				
				String[] split = inputline.split(",");
				String chromosome = split[1];
				String marker_id = split[0];
				
				
				if(!this.axiomMarkerPositions.containsKey(chromosome)) {
					this.axiomMarkerPositions.put(chromosome, new HashMap<Integer, String>());
					this.availablePositions.put(chromosome, new HashMap<String, Vector<Integer>>() );
				}
				try {
					int pos = Integer.parseInt(split[2]);
					
					this.axiomMarkerPositions.get(chromosome).put(pos, marker_id);
					
					if(!this.availablePositions.get(chromosome).containsKey(accession)) {
						this.availablePositions.get(chromosome).put(accession, new Vector<Integer>() );
					}
					
					this.availablePositions.get(chromosome).get(accession).add(pos);
					
					if(!this.axiomData.containsKey(marker_id)) {
						this.axiomData.put(marker_id, new HashMap<String, String>());
					}
					for( int i = 5; i< header.length; i++) {
						this.axiomData.get(marker_id).put(header[i], split[i]);
					}
				}catch(NumberFormatException e) {}
				
			}

			in.close();
		}
		
		
		
		
		
		for(Iterator<String> iterator = this.availablePositions.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			for(Iterator<String> iterator2 = this.availablePositions.get(chr).keySet().iterator(); iterator2.hasNext();) {
				String accession = iterator2.next();
				Collections.sort(this.availablePositions.get(chr).get(accession));
			}
			
		}
		
		
	}
	
	
	/**
	 * 
	 * Read in the table for converting the prefix of the axiom file names to the accession names that can be searched for in the vcf files
	 * 
	 * @param nameConversionTableFile
	 * 			File with the mappings. Got that from Luzie Wingen. Prefix is in column 0 and WatSeq accession name in column 5.
	 * @return
	 * @throws IOException
	 */
	private HashMap<String, String> readNameConversionTable(File nameConversionTableFile)throws IOException{
		HashMap<String, String> nameConversionTable = new HashMap<String, String>();
		BufferedReader in = new BufferedReader(new FileReader(nameConversionTableFile));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			nameConversionTable.put(split[0], split[5]);
		}

		in.close();
		return nameConversionTable;
		
	}
	
	private void readNucmerAlignments(File nucmerAlignmentFile)throws IOException{
		if(this.nucmerAlignments == null) {
			this.nucmerAlignments = new HashMap<String, Vector<int[]>>();
		}
		
		BufferedReader in = new BufferedReader(new FileReader(nucmerAlignmentFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			try {
				int[] a = {Integer.parseInt(split[0]), Integer.parseInt(split[1])};
				String chr = split[split.length-2].split("_")[0];
				chr = chr.replaceAll("chr", "");
				if( !this.nucmerAlignments.containsKey(chr)) {
					this.nucmerAlignments.put(chr, new Vector<int[]>());
				}
				this.nucmerAlignments.get(chr).add(a);
			}catch(NumberFormatException e) {} 
			
		}

		in.close();
		
		for(Iterator<String> iterator = this.nucmerAlignments.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			
			Collections.sort(this.nucmerAlignments.get(chr), new Comparator<int[]>() {public int compare(int[] a1, int[] a2) {
																		if(a1[0] < a2[0]) {return -1;} 
																		else if(a1[0] > a2[0]) {return 1;} 
																		else {return 0;}
													}});
			
		}
		
		
	}
	
	
	private void readNucmerSNPs(File nucmerSnpFile)throws IOException{
		if (this.nucmerSnps == null) {
			this.nucmerSnps = new HashMap<String, HashMap<Integer, String[]>>();
		}
		
		
		BufferedReader in = new BufferedReader(new FileReader(nucmerSnpFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\\s+");
			try {
				int pos = Integer.parseInt(split[0]);
				String chr = split[split.length-2].split("_")[0].replaceAll("chr", "");
				String[] allele = {split[1] , split[2]};//the first is the reference allele as control
				if(!this.nucmerSnps.containsKey(chr)) {
					this.nucmerSnps.put(chr, new HashMap<Integer, String[]> ());
				}
				this.nucmerSnps.get(chr).put(pos, allele);
			}catch(NumberFormatException e) {}
			
		}

		in.close();
		
	}
	
	/**
	 * 
	 * @deprecated
	 * 
	 * @param nucmer_snps
	 * @param nucmer_alignments
	 * @param vcfInput
	 * @param vcfOutput
	 * @throws IOException
	 */
	public static void mergeFiles(File nucmer_snps,File nucmer_alignments,File vcfInput,File vcfOutput)throws IOException{
		
		HashMap<String, Vector<int[]>> alignments = new HashMap<String, Vector<int[]>>();
		HashMap<String, HashMap<Integer, String>> snps = new HashMap<String, HashMap<Integer, String>>();
		
		BufferedReader in = new BufferedReader(new FileReader(nucmer_alignments));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			try {
				int[] a = {Integer.parseInt(split[0]), Integer.parseInt(split[1])};
				String chr = split[split.length-2].split("_")[0];
				
				if( !alignments.containsKey(chr)) {
					alignments.put(chr, new Vector<int[]>());
				}
				alignments.get(chr).add(a);
			}catch(NumberFormatException e) {} 
			
		}

		in.close();
		
		for(Iterator<String> iterator = alignments.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			
			Collections.sort(alignments.get(chr), new Comparator<int[]>() {public int compare(int[] a1, int[] a2) {
																		if(a1[0] < a2[0]) {return -1;} 
																		else if(a1[0] > a2[0]) {return 1;} 
																		else {return 0;}
													}});
			
		}
		
		in = new BufferedReader(new FileReader(nucmer_snps));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			try {
				int pos = Integer.parseInt(split[0]);
				String chr = split[split.length-2].split("_")[0];
				String allele = split[1] + "\t" + split[2];
				if(!snps.containsKey(chr)) {
					snps.put(chr, new HashMap<Integer, String> ());
				}
				snps.get(chr).put(pos, allele);
			}catch(NumberFormatException e) {}
			
		}

		in.close();
		
		
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(vcfOutput))));
		
		 in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream (vcfInput))));

		 int numSNPs = 0;
		 int numAligned = 0;
		 
		HashMap<String, Integer> accessionCollumns = new HashMap<String, Integer>();
		 
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			if(inputline.startsWith("#")) {
				out.write(inputline);
				
				if(inputline.startsWith("#CHROM")) {
					out.write("\tParagon");
					
					String[] split = inputline.split("\t");
					for( int i = 0; i< split.length; i++) {
						accessionCollumns.put( split[i], i);
					}
					
				}
				
				out.newLine();
			}else {
				String[] split = inputline.split("\t");
				String chr = split[0];
				int pos = Integer.parseInt(split[1]);
				String referenceAllele = split[3] ;
				String alternativeAllele = split[4];
				String allele = referenceAllele + "\t" + alternativeAllele;
				//String genotypeFieldDescription = split[8];
				
				
				if( snps.containsKey(chr) && snps.get(chr).containsKey(pos)) {
					if( allele.equalsIgnoreCase(snps.get(chr).get(pos))) {
						numSNPs++;
					}else {
						System.out.println("alleles not matching: " +pos + "\t"+ allele + "\t" + snps.get(chr).get(pos));
					}
				}else {
					
					if(snps.containsKey(chr)) {
						for(Iterator<int[]> iterator = alignments.get(chr).iterator(); iterator.hasNext();) {
							int[] a = iterator.next();
							if( pos >= a[0] && pos <= a[1]) {
								numAligned++;
							}else {
								if( a[0] > pos) {
									break;
								}
							}
							
						}
						
					}
					
				}
				
			}
			
			
		}

		in.close();
		
		out.close();
		System.out.println("number of snps: " + numSNPs);
		System.out.println("number of snps same to CS: " + numAligned);
		System.out.println(numSNPs + numAligned);
	}
	
}
