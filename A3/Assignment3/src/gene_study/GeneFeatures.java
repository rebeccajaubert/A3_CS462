package gene_study;



import java.util.LinkedHashMap;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceMixin;
import org.biojava.nbio.core.sequence.template.SequenceView;

public class GeneFeatures {
	
	private static final NucleotideCompound A = DNACompoundSet.getDNACompoundSet().getCompoundForString("A");
	private static final NucleotideCompound T = DNACompoundSet.getDNACompoundSet().getCompoundForString("T");
	private static final NucleotideCompound C = DNACompoundSet.getDNACompoundSet().getCompoundForString("C");
	private static final NucleotideCompound G = DNACompoundSet.getDNACompoundSet().getCompoundForString("G");
	private static final String[] DNA_CODONS = {"GCT", "GCC", "GCA", "GCG", "TGT", "TGC","TGG","GAT", "GAC",
		    "GAA", "GAG","TTT", "TTC","GGT", "GGC", "GGA", "GGG","CAT", "CAC","ATA", "ATT", "ATC",
		    "AAA", "AAG", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATG","AAT", "AAC","CCT", "CCC", "CCA",
		    "CCG","CAA", "CAG","CGT", "CGC", "CGA", "CGG", "AGA", "AGG","TCT", "TCC", "TCA", "TCG", "AGT", 
		    "AGC", "ACT", "ACC", "ACA", "ACG","GTT", "GTC", "GTA", "GTG","TAT", "TAC", "TAA", "TAG", "TGA" };
	
	
	
	public static LinkedHashMap[] getFrequencyTable(LinkedHashMap<String, DNASequence> dnaSequences, FeatureList features){
		
		int countAinter=0, countTinter=0, countCinter=0, countGinter=0;
		int countIntergenicNucleotides=0;
		int nbEmCodons =0;int nbStartCodons =0;int nbStopCodons =0;
		LinkedHashMap[] frequencies = new LinkedHashMap[3];
		LinkedHashMap<DNASequence, Float> startFrequencies = new LinkedHashMap<DNASequence, Float>();
		LinkedHashMap<DNASequence, Float> emissionFrequencies = new LinkedHashMap<DNASequence, Float>();	//hold the # of times a codon appeared
		LinkedHashMap<DNASequence, Float> stopFrequencies = new LinkedHashMap<DNASequence, Float>();
		
		//from string of nucleotides create dnaSequence
		for(String c : DNA_CODONS) { 
			try { 
				DNASequence seq = new DNASequence(c);
				emissionFrequencies.putIfAbsent(seq, 0f);	//initialize map of frequencies
				startFrequencies.putIfAbsent(seq, 0f);
				stopFrequencies.putIfAbsent(seq, 0f);
			} catch (CompoundNotFoundException e) {e.printStackTrace();} 
		}
		
		for (String key : dnaSequences.keySet()) {
			
			//first get the features within same contig
			int indexSpace = key.indexOf(" "); 
			final String contigString = key.substring(0, indexSpace);  //  first word of key : format "DN38.contig%d"
			FeatureList sameContigList = new FeatureList();
			for (FeatureI feature : features) {
				if(feature.seqname().equals(contigString)) {
					sameContigList.add(feature);
				}
				else if(!sameContigList.isEmpty()) { 	//we found all the corresponding features no need to keep searching
					break;
				}
			}
			
			Location prevLoc = new Location(1, 1);
			for(FeatureI feature : sameContigList) {
			
				//previous intergenic region
				SequenceView<NucleotideCompound> intergenicSequence = SequenceMixin.createSubSequence(dnaSequences.get(key), prevLoc.bioEnd(), feature.location().bioStart());
				for(NucleotideCompound nucleotideCompound : intergenicSequence) {	//count all nucleotides
					if(nucleotideCompound.equals(A) ) countAinter++;
					if(nucleotideCompound.equals(T) ) countTinter++;
					if(nucleotideCompound.equals(C) ) countCinter++;
					if(nucleotideCompound.equals(G) ) countGinter++;
					countIntergenicNucleotides++;
				}
				
				//genic region
				SequenceView<NucleotideCompound> genicSequence = SequenceMixin.createSubSequence(dnaSequences.get(key), feature.location().bioStart(), feature.location().bioEnd());
				int codonCounter = 0; //if =3 then start new codon
				String codonString ="";	//hold the current codon
				int counterHelperStartStop =0;
				for(NucleotideCompound nucleotideCompound : genicSequence) {
					
					if(codonCounter<2) {
						codonString += nucleotideCompound.toString();
						codonCounter++;
						counterHelperStartStop++;
					}
					else {
						codonString += nucleotideCompound.toString();
						try { 
							DNASequence codonSeq = new DNASequence(codonString);
							if(counterHelperStartStop == 2) { 	//start codon
								startFrequencies.put(codonSeq,startFrequencies.get(codonSeq)+1);
								nbStartCodons++;
							}
							else if(counterHelperStartStop == genicSequence.getLength()-1) {	//stop codon
								stopFrequencies.put(codonSeq,stopFrequencies.get(codonSeq)+1);
								nbStopCodons++;
							}
							else {
								emissionFrequencies.put(codonSeq, emissionFrequencies.get(codonSeq) + 1); //increment
								nbEmCodons++;
							}
						}catch (CompoundNotFoundException e) {e.printStackTrace();} 
						
						codonCounter=0; codonString=""; //new codon
						counterHelperStartStop++;
					}
				}
				prevLoc = feature.location();
			}
		}
		
		final int totalEmCodons = nbEmCodons; final int totalStartCodons = nbStartCodons; final int totalStopCodons = nbStopCodons;
		//compute frequency of each codon
		emissionFrequencies.replaceAll( (k,v)-> v = v/(float)totalEmCodons );
		startFrequencies.replaceAll( (k,v)-> v = v/(float)totalStartCodons );
		stopFrequencies.replaceAll( (k,v)-> v = v/(float)totalStopCodons );
		
		float freqA = (float)countAinter/(float)countIntergenicNucleotides;
		float freqT = (float)countTinter/(float)countIntergenicNucleotides;
		float freqC = (float)countCinter/(float)countIntergenicNucleotides;
		float freqG = (float)countGinter/(float)countIntergenicNucleotides;
		//add freq of A,T,C,G  to hashmap
		try {
			DNASequence seqA = new DNASequence("A");
			emissionFrequencies.putIfAbsent(seqA, freqA);
			DNASequence seqT = new DNASequence("T");
			emissionFrequencies.putIfAbsent(seqT, freqT);
			DNASequence seqC = new DNASequence("C");
			emissionFrequencies.putIfAbsent(seqC, freqC);
			DNASequence seqG = new DNASequence("G");
			emissionFrequencies.putIfAbsent(seqG, freqG);
		}catch (CompoundNotFoundException e) {e.printStackTrace();} 
		
		frequencies[0] = startFrequencies;
		frequencies[1] = emissionFrequencies;
		frequencies[2] = stopFrequencies;

		return  frequencies;
	}
	
	
	
	
	
	//Compute both intergenic and genic regions length in the same method for an efficiency matter
	public static Pair<Float> getAverageLengthRegions(FeatureList features) {
		int sumIntergenicLengths =0;
		int sumGenicLengths =0;
		
		//iterate through genes (counter position restart at each gene)
		for (int k=0; k< features.size()-1; k++){
			//first gene sequence with a new seqname
			Feature prevFeature = (Feature) features.get(k);
			Location prevFeatureLocation = prevFeature.location();
			//initialize the var holding the summations of every intergenic/genic regions' length 
			sumIntergenicLengths += prevFeatureLocation.bioStart();				//first intergenic region if start location is not 0
			sumGenicLengths += prevFeatureLocation.bioEnd()- prevFeatureLocation.bioStart(); //first genic length 
			
			k++;
			Feature feature =  (Feature) features.get(k);
			while(prevFeature.seqname().equals(feature.seqname())){ //same sequence
				Location featureLocation = feature.location();
				
				//genic length
				int geneLength = featureLocation.bioEnd() - featureLocation.bioStart();
				sumGenicLengths += geneLength;
				
				//intergenic length
				int interLength = featureLocation.bioStart() - prevFeatureLocation.bioEnd();	//length of the intergenic region between genes i and i-1
				if(interLength<0) { //in case genes overlap then ignore
					prevFeatureLocation = featureLocation;
					k++;
					feature =  (Feature) features.get(k);
					continue;
				}
				sumIntergenicLengths += interLength;
				prevFeatureLocation = featureLocation;
				
				k++;
				feature =  (Feature) features.get(k);
			}
		}
		//compute average lengths
		float intergenicAvg = sumIntergenicLengths / features.size();
		float genicAvg = sumGenicLengths / features.size();
		Pair<Float> averages = new Pair<Float>(intergenicAvg, genicAvg); 
		return averages;
	}

}
