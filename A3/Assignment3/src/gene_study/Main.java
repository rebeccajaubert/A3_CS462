package gene_study;

import gene_study.GeneFeatures;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.LinkedHashMap;
import java.util.function.Predicate;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.structure.contact.Pair;
import org.forester.applications.get_distances;

public class Main {
	public static void main(String[] args) {

		System.out.println("Gene study");

		try  
		{  
			FeatureList featuresList = GFF3Reader.read("./files/Vibrio_cholerae.GFC_11.37.gff3");	//read file to get all genes
			FeatureList featuresListQ1D = GFF3Reader.read("./files/Vibrio_vulnificus.ASM74310v1.37.gff3");
			
			//we consider only the Coding Sequences and positive strands
			featuresList.removeIf( new Predicate<FeatureI>() {
				public boolean test(FeatureI g) {
					return g.location().isNegative()  ||  !g.type().equals("CDS");
				}
			} );
			
			featuresListQ1D.removeIf( new Predicate<FeatureI>() {
				public boolean test(FeatureI g) {
					return g.location().isNegative()  ||  !g.type().equals("CDS");
				}
			} );
			
			Pair<Float> lengths = GeneFeatures.getAverageLengthRegions(featuresList);
//			System.out.println("Average length of the intergenic regions : " + lengths.getFirst() );
//			System.out.println("Average length of the genic regions : " + lengths.getSecond() );

			File Q1Cfile = new File("files/Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa");
			File fastaFile = new File("files/Vibrio_cholerae.GFC_11.dna.toplevel.fa");
			LinkedHashMap<String, DNASequence> dnaSequences = FastaReaderHelper.readFastaDNASequence(fastaFile,true);
			LinkedHashMap<String, DNASequence> dnaSequencesQ1C = FastaReaderHelper.readFastaDNASequence(Q1Cfile,true);
			
			LinkedHashMap[] frequencies = GeneFeatures.getFrequencyTable(dnaSequences, featuresList);

			//create file with all the information found in a and probabilities
			File infoFile = new File("./files/information.txt");
			FileOutputStream fos = new FileOutputStream(infoFile);
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
			bw.write("length averages : Intergenic , Genic"); bw.newLine();
			bw.write(lengths.getFirst().toString());bw.newLine(); //first line is avg length intergenic region
			bw.write(lengths.getSecond().toString());bw.newLine(); // second line is avg length genic region
			
			for(int i=0; i<frequencies.length;i++) {
				for(Object sequence : frequencies[i].keySet()) {
					try {
						bw.write(sequence.toString() + " " + frequencies[i].get(sequence)); bw.newLine();
					} catch (IOException e) {
						e.printStackTrace();
					}  
				}
				bw.newLine();
			}

			//write probabilities
			bw.write("Probabilities: Intergenic/Intergenic , Inter/Start, Start/Gene, Gene/Gene, Gene/Stop, Stop/Inter "); bw.newLine();
			
			
			float interInterTrans = (lengths.getFirst()-1)/lengths.getFirst();
			float interStartTrans = 1/lengths.getFirst();
			float geneGeneTrans = (lengths.getSecond()/3-1)/(lengths.getSecond()/3);
			float geneStopTrans = 3/lengths.getSecond();
			
			float[] transitions = {interInterTrans,interStartTrans,1f,geneGeneTrans,geneStopTrans,1f};
			for(int i=0; i<transitions.length; i++) {
				bw.write(Float.toString(transitions[i]));bw.newLine();
			}

			bw.close();
			
			File genesFile = new File("./files/genes.txt");
			
			//to run question question 1C change "dnaSequences" to "dnaSequencesQ1C"
			//GenePrediction.Viterbi(fastaFile, infoFile, dnaSequencesQ1C, frequencies[0], frequencies[1], frequencies[2], transitions, genesFile);
			
			
			
			File genesFile1C = new File("files/1Cgenes.txt");
			AccurcyPredictions.getAccuracy(featuresListQ1D, genesFile1C);
			
		}
		catch (IOException ex) {
			ex.printStackTrace();
		}
	}

}
