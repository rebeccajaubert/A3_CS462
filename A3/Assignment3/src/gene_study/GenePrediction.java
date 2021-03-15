package gene_study;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.LinkedHashMap;


import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;

public class GenePrediction {
	
	private static String[] states = { "I", "Start", "G", "Stop"};
	private static float[] initial_probabilities = {1,0,0,0};
	private static float[][] transition_probabilities = { {0,0},{0,0},{0,0},{0,0} }; //4 states : Intergenic, start, gene, stop //first entry = state to same state , 2nd entry = state to next //no need for 4 entries cause would be zeros
	
	
	private static class Node {
		public float probGenic;	
		public float probStart;
		public float probIntergenic;	
		public float probStop;
		
		public DNASequence nuclDnaSequence;
		public DNASequence codonDnaSequence;
		
		public Node( float intergenicProb, float startProb, float genicProb, float stopProb,DNASequence currentNucl ,DNASequence currentCodon) {
			this.probGenic = genicProb;
			this.probStart = startProb;
			this.probIntergenic = intergenicProb;
			this.probStop = stopProb;
			this.nuclDnaSequence = currentNucl;
			this.codonDnaSequence = currentCodon;
		}
		
		@Override
		public String toString() {
			return " " + probIntergenic + " " + probStart + " " + probGenic + " " + probStop;
		}
	}
	
	private static float log(float a) {
		if(a==0) return -9999;
		else return (float) Math.log(a);
	}
	
	public static void Viterbi(File fasta, File info, LinkedHashMap<String, DNASequence> dnaSequences,
			LinkedHashMap<DNASequence, Float>   startFrequencies,LinkedHashMap<DNASequence, Float> emissionFrequencies,
			LinkedHashMap <DNASequence, Float> stopFrequencies ,float[] transitions, File genesFile) {
		
		/*
		 * The file info in the arguments contains all the other arguments but it was faster to reuse directly my structures
		 * obtained in q1a by running methods in GeneFeatures instead of recreating them by reading the file
		 */
		
		
		//output genes in genesFile
		BufferedWriter bw = null ;
		try {
			FileOutputStream fos2 = new FileOutputStream(genesFile);
			 bw = new BufferedWriter(new OutputStreamWriter(fos2));
		}catch (IOException e) {
			e.printStackTrace();
		} 
		
		
		transition_probabilities[0][0] = transitions[0]-0.1f; transition_probabilities[0][1] = transitions[1]+0.1f; //it works without modif but better accuracy
		transition_probabilities[1][1] = transitions[2];
		transition_probabilities[2][0] = transitions[3]; transition_probabilities[2][1] = transitions[4];
		transition_probabilities[3][1] = transitions[5];
		
		for (String contig : dnaSequences.keySet()) {
			
			DNASequence currentSequence = dnaSequences.get(contig);
			String sequenceAsString = currentSequence.getSequenceAsString();
			
			//initialization matrix
			Node[] V = new Node[sequenceAsString.length()]; //Node will act as an array 

			//iterate through nucleotides
			for(int j=2; j< sequenceAsString.length() ;j++) {

				String curNuclString = sequenceAsString.charAt(j) + "";
				String prevNuclString = sequenceAsString.charAt(j-1) + "";
				String prevPrevNuclString = sequenceAsString.charAt(j-2) + "";
				String curCodon = curNuclString+prevNuclString+prevPrevNuclString;
				DNASequence currentNucl = null; //DNASequence nextNucl = null; DNASequence nextNextNucl = null;
				DNASequence currentCodon = null;
				
				try { //easier to deal with a nucleotide as a dnasequence 
					currentNucl = new DNASequence(curNuclString);
					currentCodon = new DNASequence(curCodon);
					
					if(j==2) {	//initialize
						float genicProb =  Math.max((log(initial_probabilities[2]) + log(emissionFrequencies.get(currentCodon))),
									(log(initial_probabilities[1]) + log(emissionFrequencies.get(currentCodon))) );

						float interProb =  Math.max((log(initial_probabilities[0]) + log(emissionFrequencies.get(currentNucl))),
										(log(initial_probabilities[3]) + log(emissionFrequencies.get(currentCodon))));
						float startProb =  (log(initial_probabilities[1]) + log(startFrequencies.get(currentCodon)));
						float stopProb =  (log(initial_probabilities[3]) + log(stopFrequencies.get(currentCodon)));
						
						Node entryV = new Node(interProb, startProb, genicProb, stopProb, currentNucl, currentCodon); //first entry V
						V[j] = entryV; V[j-1] = entryV; V[j-2] = entryV; 
						
						//System.out.println(entryV.toString());
						continue;
					}

					float interProb = Math.max( (V[j-1].probIntergenic) + log(emissionFrequencies.get(currentNucl)) +  log(transition_probabilities[0][0]) ,
								(V[j-1].probStop) + log(emissionFrequencies.get(currentCodon)) + log(transition_probabilities[3][1]) );
					float genicProb =  Math.max( (V[j-3].probStart) + log(emissionFrequencies.get(currentCodon)) + log(transition_probabilities[1][1]) ,
									(V[j-3].probGenic) + log(emissionFrequencies.get(currentCodon)) + log(transition_probabilities[2][0]) ) ;
					float startProb =  ((V[j-1].probIntergenic) + log(startFrequencies.get(currentCodon)) + log(transition_probabilities[0][1]));
					float stopProb =  ((V[j-3].probGenic) + log(stopFrequencies.get(currentCodon)) + log(transition_probabilities[2][1]));
					
					Node entryV = new Node(interProb, startProb, genicProb, stopProb, currentNucl, currentCodon);
										
					V[j] = entryV;
					//System.out.println(entryV.toString());
					
				}
				catch (CompoundNotFoundException e) {	e.printStackTrace();}

				
			}//finished iterating through nucleotide sequence
			
			//traceback
			String path = ""; int indexState=-1; 
			int stopIndex=0;
			String[] contigInfo = contig.split(" "); 
			
			//init traceback
			//find where max comes from
			float max = Math.max(  Math.max(V[sequenceAsString.length()-1].probStart,V[sequenceAsString.length()-1].probIntergenic ) ,
								Math.max(V[sequenceAsString.length()-1].probGenic, V[sequenceAsString.length()-1].probStop));
			if (max == V[sequenceAsString.length()-1].probIntergenic) 
				indexState=0;
			if( max==V[sequenceAsString.length()-1].probStart)
				indexState=1;
			if( max == V[sequenceAsString.length()-1].probGenic)
				indexState=2;
			if (max == V[sequenceAsString.length()-1].probStop)
				indexState=3;
			
			String curState = states[indexState];
			
			//find how we got in this state
			for(int j=sequenceAsString.length()-2 ;j>1; j--) {
				path+=curState;
				if(curState == states[0]) { //inter
					if( V[j].probIntergenic == V[j-1].probIntergenic+log(emissionFrequencies.get(V[j].nuclDnaSequence))+log(transition_probabilities[0][0]) )  //from Inter
						curState=states[0]; 
					else { 
						curState = states[3];
						
					}	//from stop
				}
				else if (curState == states[1]) { //start
					curState=states[0]; //has to come from inter
					try {
						bw.write(contigInfo[0] + "		" + j + " " + stopIndex );	//write to file
						bw.newLine();
					} catch (IOException e) {
						e.printStackTrace();
					} 
					
					
				}
				else if (curState == states[2]) { //gene
					if( V[j].probGenic == V[j-3].probGenic+log(emissionFrequencies.get(V[j].codonDnaSequence))+log(transition_probabilities[2][0]) )   //from gene
						curState=states[2];
					else 
						curState = states[1];	//from Start
				}
				else {	//stop
					curState=states[2]; //has to come from gene
					stopIndex=j;
				}
			}
		}
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}



}

