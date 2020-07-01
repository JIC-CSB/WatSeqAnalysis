

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import support.BioSequence;
import support.Kmer;
import support.commandLineInterface.CLI;
import support.commandLineInterface.CLIParseException;

public class AxiomKmerSet {

	
	HashMap<BitSet, String> kmers;
	
	HashMap<Integer, String> accession_to_number;
	HashMap<String, HashSet<Integer>> allele_to_accessions;
	
	public final int kmer_size = 31;
	
	public static void main(String[] args) {
		
		CLI cli = new CLI();
		
		cli.parseOptions(args);
		
		String help = ("java -jar CheckAxiomData.jar -m axiom_manifest.tsv -a axiom_data.tsv -i kmer_dump -o output_file.txt" + "\n"  
					+ "-m\tthe axiom manifest file as tsv. The ID is considered to be in the first column, the sequence in the fifth. " + "\n"		
					+ "-a\tthe axiom data. First row contains the accession names. Columns with genotypes for accessions start from (including) 12th column." + "\n"
					+ "-i\tthe kmer dump as tsv. kmers are in the first column. k=31" + "\n"
					+ "-o\toutput prefix" + "\n"
				);
		
		try {
			
			if( !cli.hasOption("m") || !cli.hasOption("a") || !cli.hasOption( "i") || !cli.hasOption("o")){
				throw new CLIParseException(" missing parameters!");
				
			}
				
			
			File axiomManifest = new File(cli.getArg("m"));
			File axiomWatkins = new File(cli.getArg("a"));
			File kmerDump = new File(cli.getArg("i"));
			File outputFileResult = new File(cli.getArg("o") + ".result.txt");
			File outputFileFilteredDump = new File(cli.getArg("o") + ".filtered_dump.txt");
			File inSilicoGenotype = new File(cli.getArg("o") + ".genotype.txt");
			
			AxiomKmerSet axiom = new AxiomKmerSet(axiomManifest,axiomWatkins);
			
			HashMap<BitSet, Integer> kmerTable = axiom.readKmerSet(kmerDump);
			
			axiom.writeFilteredKmerSet(kmerTable, outputFileFilteredDump);
			
			HashMap<String, HashMap<String, int[]>> genotypes = axiom.in_silico_genotyping(kmerTable, inSilicoGenotype);
			HashMap<String ,HashMap<String,String>> axiomGenotypes = axiom.getAxiomGenotypes();
			
			
			axiom.getGenotypeComparison(genotypes, axiomGenotypes, outputFileResult);
			
			
			
			
			
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println(help);
		}
		catch( CLIParseException e) {
			e.printStackTrace();
			System.out.println(help);
		}
		
		
	}
	
	
	public static void mainTest(String[] args) {
		try {
			File axiomManifest = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/SimonGriffiths/WatSeq/axiom/35k_probe_set_IWGSCv1-1.txt");
			File axiomWatkins = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/SimonGriffiths/WatSeq/axiom/Watkins_829acc_axiom35k_gts_hmp.txt");
			
			
			File kmerDump = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/SimonGriffiths/WatSeq/axiom/kmerData/WATDE0703/accession.filtered_dump.txt");
			File outputFile = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/SimonGriffiths/WatSeq/axiom/kmerData/test.txt");
			
			
			AxiomKmerSet axiom = new AxiomKmerSet(axiomManifest,axiomWatkins);
			
			HashMap<BitSet, Integer> kmerTable = axiom.readKmerSet(kmerDump);
			System.out.println("kmers read");
			
			//HashMap<String, Integer> accessionCount = axiom.checkDumpIncludeNegative(kmerTable);
			
			
			
			
			//HashMap<String, Integer> accessionCount = axiom.checkDump(kmerDump);
			//axiom.writeAccessionCount(accessionCount, outputFile);
			
			HashMap<String, HashMap<String, int[]>> genotypes = axiom.in_silico_genotyping(kmerTable, outputFile);
			HashMap<String ,HashMap<String,String>> axiomGenotypes = axiom.getAxiomGenotypes();
			
			
			axiom.getGenotypeComparison(genotypes, axiomGenotypes, outputFile);
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	
	
	
	
	public AxiomKmerSet(File axiomManifest, File axiomWatkins)throws IOException{
		
		//initialize data structures
		this.kmers = new HashMap<BitSet, String>();
		this.accession_to_number = new HashMap< Integer, String>();
		this.allele_to_accessions = new HashMap<String, HashSet<Integer>>();
		
		
		//read in data
		readAxiomManifest(axiomManifest);
		readAxiomWatkinsData(axiomWatkins);
	}
	
	
	private void readAxiomManifest(File axiomManifest)throws IOException{
		
		HashSet<BitSet> remove = new HashSet<BitSet>();
		
		BufferedReader in = new BufferedReader(new FileReader(axiomManifest));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			String id = split[0];
			String sequence = split[4];
			//TGTTCA[A/G]GCGTTA
			String s1 = sequence.split("\\[")[0] ;
			String s2 = sequence.split("\\]")[1];
			String snp = sequence.split("\\[")[1].split("\\]")[0];
			
			int k1 = Math.min(s1.length(),kmer_size-1 );
			int k2 = Math.min(s2.length(),kmer_size-1 );
			
			String sequence1 = s1.substring(s1.length()-k1) + snp.substring(0,1) + s2.substring(0,k2);
			String sequence2 = s1.substring(s1.length()-k1) + snp.substring(2,3) + s2.substring(0,k2);
			
			
			for( int i = 0; i< sequence1.length()- kmer_size; i++) {
				BitSet kmer = Kmer.convert(sequence1.substring(i, i+kmer_size));
				
				if( this.kmers.containsKey(kmer)) {
					//System.out.println("Kmer existing\t" + kmers.get(kmer) + "\t" +id + "\t" + snp.substring(0,1) );
					remove.add(kmer);
				}
				
				this.kmers.put(kmer, id + "\t" + snp.substring(0,1));
			}
			
			for( int i = 0; i< sequence2.length()-kmer_size; i++) {
				BitSet kmer = Kmer.convert(sequence2.substring(i, i+kmer_size));
				
				if( this.kmers.containsKey(kmer)) {
					//System.out.println("Kmer existing\t" + kmers.get(kmer) + "\t" +id + "\t" + snp.substring(2,3) );
					remove.add(kmer);
				}
				
				this.kmers.put(kmer, id + "\t" + snp.substring(2,3));
			}
			
		}
		
		in.close();
		System.out.println("Kmers recorded " + this.kmers.size());
		
		for( Iterator<BitSet> iterator = remove.iterator(); iterator.hasNext();) {
			this.kmers.remove(iterator.next());
		}
		System.out.println("Kmers after filtering " + this.kmers.size());
		
	}
	
	
	private void readAxiomWatkinsData(File axiomWatkins)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(axiomWatkins));
		String[] header = in.readLine().split("\t");
		
		for( int i = 11; i< header.length; i++ ) {
			this.accession_to_number.put( i,header[i]);
		}
		
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			String id = split[0];
			for( int i = 11; i< split.length; i++) {
				String s = split[i];
				if(s.equalsIgnoreCase("NN")) {
					continue;
				}
				String id1 = id + "\t" + s.substring(0,1);
				String id2 = id + "\t" + s.substring(1,2);
				
				if( !allele_to_accessions.containsKey(id1)) {allele_to_accessions.put(id1, new HashSet<Integer>());}
				if( !allele_to_accessions.containsKey(id2)) {allele_to_accessions.put(id2, new HashSet<Integer>());}
				allele_to_accessions.get(id1).add(i);
				allele_to_accessions.get(id2).add(i);
			}
			
			
		}

		in.close();
		
		//System.out.println(allele_to_accessions.size());
		
	}
	/**
	 * @deprecated
	 * 
	 * @param dumpFile
	 * @return
	 * @throws IOException
	 */
	public HashMap<String, Integer> checkDump(File dumpFile)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(dumpFile));
		HashMap<String, Integer> accessionCount = new HashMap<String, Integer>();
		
		int counta = 0;
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			
			BitSet kmer1 = Kmer.convert(split[0]);
			BitSet kmer2 = Kmer.convert(new BioSequence("",split[0]).getReverseComplementarySequence());
			if( this.kmers.containsKey(kmer1)) {
				String allele = kmers.get(kmer1);
				if( ! allele_to_accessions.containsKey(allele)) {
					continue;
				}
				
				HashSet<Integer> accessions = allele_to_accessions.get(allele);
				for(Iterator<Integer> iterator = accessions.iterator(); iterator.hasNext();) {
					int i = iterator.next();
					String accession = accession_to_number.get(i);
					int num=0;
					if( accessionCount.containsKey(accession)) {
						num=accessionCount.get(accession);
					}
					num++;
					accessionCount.put(accession, num);
				}
			}
			if( this.kmers.containsKey(kmer2)) {
				String allele = kmers.get(kmer2);
				if( ! allele_to_accessions.containsKey(allele)) {
					continue;
				}
				HashSet<Integer> accessions = allele_to_accessions.get(allele);
				for(Iterator<Integer> iterator = accessions.iterator(); iterator.hasNext();) {
					int i = iterator.next();
					String accession = accession_to_number.get(i);
					int num=0;
					if( accessionCount.containsKey(accession)) {
						num=accessionCount.get(accession);
					}
					num++;
					accessionCount.put(accession, num);
				}
			}
			counta ++;
			if(counta%1000000==0) {
				System.out.print(".");
				if(counta%100000000==0) {
					System.out.println();
				}
			}
		}

		in.close();
		return accessionCount;
	}
	
	
	public HashMap<String, Integer> checkDumpIncludeNegative(File kmerDump)throws IOException{
		
		HashMap<BitSet, Integer> kmerTable = this.readKmerSet(kmerDump);
		return this.checkDumpIncludeNegative(kmerTable);
		
	}
	
	public HashMap<String, Integer> checkDumpIncludeNegative(HashMap<BitSet, Integer> kmerTable){
		HashMap<String, Integer> accessionCount = new HashMap<String, Integer>();
		
		
		for(Iterator<Integer> iterator = this.accession_to_number.keySet().iterator(); iterator.hasNext();) {
			int num = iterator.next();
			accessionCount.put(this.accession_to_number.get(num),0);
		}
		
		
		
		for (Iterator<BitSet> iteratorKmers = kmerTable.keySet().iterator(); iteratorKmers.hasNext();) {
			
			
			BitSet kmer1 = iteratorKmers.next();
			
			String allele = null;
			if(kmers.containsKey(kmer1)) {
				allele = kmers.get(kmer1);
			}else {
				BitSet kmer2 = Kmer.complement(kmer1, 31);
				if( kmers.containsKey(kmer2)) {
					allele = kmers.get(kmer2);
				}
			}
			
			if( allele != null) {
				HashSet<Integer> accessions = allele_to_accessions.get(allele);
				for(Iterator<Integer> iterator = this.accession_to_number.keySet().iterator(); iterator.hasNext();) {
					int accession = iterator.next();
					String id = this.accession_to_number.get(accession);
					
					int num = accessionCount.get(id);
					//System.out.println(id + "\t" + accession +"\t" + num + "\t" + allele);
					if(accessions!=null && accessions.contains(accession)) {
						num++;
					}else {
						num--;
					}
					accessionCount.put(id, num);
				}
				
			}
			
		}

		
		
		String topAccession = "";
		int topCount = 0;
		for(Iterator<String> iterator = accessionCount.keySet().iterator(); iterator.hasNext();) {
			String id = iterator.next();
			int num = accessionCount.get(id);
			if(num>topCount) {
				topCount = num;
				topAccession = id;
			}
		}
		
		System.out.println("Top scoring accession is " + topAccession + " with a count of " + topCount);
		return accessionCount;
	}
	
	
	
	public void writeAccessionCount(HashMap<String, Integer> accessionCount, File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Iterator<String> iterator = accessionCount.keySet().iterator(); iterator.hasNext();) {
			String accession = iterator.next();
			out.write(accession + "\t" + accessionCount.get(accession).intValue());
			out.newLine();
		}
		
		out.close();
		
		
	}


	public void writeFilteredKmerSet(HashMap<BitSet, Integer> h, File outputFile)throws IOException {
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Iterator<BitSet> iterator = h.keySet().iterator(); iterator.hasNext();) {
			BitSet kmer = iterator.next();
			String s =  Kmer.convert(kmer, kmer_size);
			int count = h.get(kmer);
			out.write(s + "\t" + count);
			out.newLine();
		}
		
		out.close();
	}

	public void writeFilteredKmerSet(File inputDump, File outputFile)throws IOException {
		
		HashMap<BitSet, Integer> h = this.readKmerSet(inputDump);
		
		this.writeFilteredKmerSet(h, outputFile);
		
	}

	
	public HashMap<BitSet, Integer> readKmerSet(File inputDump)throws IOException {
		
		HashMap<BitSet, Integer> h = new HashMap<BitSet, Integer> ();
		
		BufferedReader in = new BufferedReader(new FileReader(inputDump));
		
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			BitSet kmer1 = Kmer.convert(split[0]);
			
			
			if( this.kmers.containsKey(kmer1)) {
				
				int count = Integer.parseInt(split[1]);
				h.put(kmer1, count);
				
			}else {
				BitSet kmer2 = Kmer.convert(new BioSequence("",split[0]).getReverseComplementarySequence());
				if(this.kmers.containsKey(kmer2)) {
					int count = Integer.parseInt(split[1]);
					h.put(kmer2, count);
				}
				
			}
			
		}
		
		in.close();
		
		return h;
		
	}

	
	
	public HashMap<String, HashMap<String, int[]>> in_silico_genotyping( HashMap<BitSet, Integer> h , File outputFile)throws IOException{
		
		HashMap<String, HashMap<String, int[]>> genotype = new HashMap<String, HashMap<String, int[]>>(); // HashMap<axiom_id, HashMap<allele, {num_kmers, count_kmers }>>
		
		for(Iterator<BitSet> iterator = h.keySet().iterator(); iterator.hasNext(); ) {
			BitSet kmer = iterator.next();
			String id = kmers.get(kmer);
			if(id==null) {
				id = kmers.get(Kmer.complement(kmer, kmer_size)); 
			}
			
			String[] split = id.split("\t"); //the id is the axiom_id + \t + allele
			String axiom_id = split[0];
			String allele = split[1];
			int num = 0;
			int count = 0;
			if(!genotype.containsKey(axiom_id)) {
				genotype.put(axiom_id, new HashMap<String, int[]>() );
			}
			
			if(genotype.get(axiom_id).containsKey(allele)) {
				num = genotype.get(axiom_id).get(allele)[0];
				count = genotype.get(axiom_id).get(allele)[1];
			}
			num++;
			count = count + h.get(kmer);
			int[] a = {num, count};
			genotype.get(axiom_id).put(allele, a);
		}
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("axiom_id\tgenotype\tallele1\tallele2\tnum_allele1\tnum_allele2\tlog10_ratio\tcount_allele1\tcount_allele2");
		out.newLine();
		
		for(Iterator<String> iterator1 = genotype.keySet().iterator(); iterator1.hasNext();) {
			String axiom_id = iterator1.next();
			Vector<String> v = new Vector<String>();
			v.addAll(genotype.get(axiom_id).keySet());
			Collections.sort(v); //alphabetically sorted alleles!
			int n1 = genotype.get(axiom_id).get(v.get(0))[0];
			int n2 = 0;
			int c1 = genotype.get(axiom_id).get(v.get(0))[1];
			int c2 = 0;
			double d = 1000;
			if( v.size()>1) {
				n2 = genotype.get(axiom_id).get(v.get(1))[0];
				c2 = genotype.get(axiom_id).get(v.get(1))[1];
				d = (double)n1 / (double)n2;
				d = Math.log10(d);
				out.write(axiom_id  +"\t"+ v.get(0)+v.get(1) + "\t" + v.get(0) +"\t" +v.get(1) +"\t" + n1 + "\t" + n2 +"\t" + d +"\t"+c1 + "\t" + c2);
			}else {
				out.write(axiom_id  +"\t"+ v.get(0)+v.get(0)+ "\t" + v.get(0) +"\t" +"NA" +"\t" + n1 + "\t" + n2 +"\t" + d+"\t"+c1 + "\t" + c2);
			}
			
			out.newLine();
			
		}
		
		out.close();
		return genotype;
	}
	
	
	
	public void getGenotypeComparison(HashMap<String, HashMap<String, int[]>> genotypes, HashMap<String, HashMap<String,String>> axiomGenotypes, File outputFile)throws IOException{
		
		int scoreHet = 1;
		int scoreMatch = 2;
		int scoreMisMatch = -2;
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for( Iterator<Integer> iteratorAccessions = this.accession_to_number.keySet().iterator(); iteratorAccessions.hasNext();) {
			int accessionnumber = iteratorAccessions.next();
			String accession = this.accession_to_number.get(accessionnumber);
			
			int score = 0;
			
			for( Iterator<String> iterator = genotypes.keySet().iterator(); iterator.hasNext(); ) {
				String axiomID = iterator.next();
				String genotype = "";
				if( genotypes.get(axiomID).size() ==1){
					genotype = genotypes.get(axiomID).keySet().iterator().next() + genotypes.get(axiomID).keySet().iterator().next();
				}else {
					Iterator<String> iterator2 = genotypes.get(axiomID).keySet().iterator();
					String allele1 = iterator2.next();
					String allele2 = iterator2.next();
					genotype = allele1 + allele2;
					if( allele1.compareTo(allele2)>0) {
						genotype = allele2+allele1;
					}
				}
				
				String axiomGenotype = axiomGenotypes.get(axiomID).get(accession);
				if(axiomGenotype == null) {
					continue;
				}
				//System.out.println(genotype + "\t" + axiomGenotype);
				if( genotype.equalsIgnoreCase(axiomGenotype)) {
					score = score + scoreMatch;
				}else if(axiomGenotype.substring(0, 1).equalsIgnoreCase(axiomGenotype.substring(1, 2)) &&
						(genotype.contains( axiomGenotype.substring(0,1)) || genotype.contains(axiomGenotype.substring(1, 2)))) {
					score = score + scoreHet;
					
				}else {
					score = score + scoreMisMatch;
				}
				
			}
			out.write(accession + "\t" + score);
			out.newLine();
		}
		
		out.close();
		
		
		
	}
	
	
	/**
	 * 
	 * re-calculate the genotypes from the genotyping of the watkins collection with the axiom35k array.
	 * 
	 * @return
	 * 		A HashMap, the keys are the axiom ids. The values are again HashMaps with Watkins accession as key and the genotype as value.
	 */
	public HashMap<String, HashMap<String, String>> getAxiomGenotypes(){
		HashMap<String, HashMap<String, String>> axiomGenotypes = new HashMap<String, HashMap<String, String>>();// axiom_id, <HashMap<accession,genotype>
		
		for(Iterator<String> iterator = this.allele_to_accessions.keySet().iterator(); iterator.hasNext();) {
			String allele = iterator.next();
			String axiom_id = allele.split("\t")[0];
			String snp = allele.split("\t")[1];
			
			if( ! axiomGenotypes.containsKey(axiom_id)) {
				axiomGenotypes.put(axiom_id, new HashMap<String, String>() );
			}
			
			HashSet<Integer> accessions = allele_to_accessions.get(allele);
			for(Iterator<Integer> iterator2 = accessions.iterator(); iterator2.hasNext();) {
				int accession = iterator2.next();
				String id = accession_to_number.get(accession);
				
				if(axiomGenotypes.get(axiom_id).containsKey(id)) {
					String snp2 = axiomGenotypes.get(axiom_id).get(id);
					String genotype = snp + snp2 ;
					if( snp2.compareTo(snp)<0 ) {
						genotype = snp2 + snp;
					}
					axiomGenotypes.get(axiom_id).put(id, genotype);
				}else {
					axiomGenotypes.get(axiom_id).put(id, snp);
				}
				
			}
		}
		
		for(Iterator<String> iterator1 = axiomGenotypes.keySet().iterator(); iterator1.hasNext();) {
			String axiom_id = iterator1.next();
			for(Iterator<String> iterator2 = axiomGenotypes.get(axiom_id).keySet().iterator(); iterator2.hasNext();) {
				String id = iterator2.next();
				String genotype = axiomGenotypes.get(axiom_id).get(id);
				if(genotype.length() ==1) {
					genotype = genotype + genotype;
					 axiomGenotypes.get(axiom_id).put(id, genotype);
				}
				//System.out.println(axiom_id + "\t" + id + "\t" + genotype);
				
			
			}
					
		}
		
		
		return axiomGenotypes;
	}

}


