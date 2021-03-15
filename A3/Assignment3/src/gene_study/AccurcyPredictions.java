package gene_study;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;

public class AccurcyPredictions {
	
	public static void getAccuracy(FeatureList featuresToCompareTo, File genesFile) {
		
		int perfectMatches=0;
		int startMatches=0;
		int endMatches=0;
		int notSameEndpointsMatches=0;
		
		int NbGenesInFile =0; boolean firstIteration=true;
		
		try {
			FileReader fr=new FileReader(genesFile);
			//StringBuffer sb=new StringBuffer();
			
			for(FeatureI feature : featuresToCompareTo) {
				BufferedReader br=new BufferedReader(fr); 
				String line;  
				while((line=br.readLine())!=null)  
				{  
					if(firstIteration) NbGenesInFile++;
					String[] contigInfo = line.split("\\s+"); // contig start end
					String contigN = contigInfo[0];
					String start = contigInfo[1]; String end = contigInfo[2];

					if(contigN.equals( feature.seqname() )){
							if(start.equals( Integer.toString(feature.location().bioStart())) &&  end.equals( Integer.toString(feature.location().bioEnd()))){
								perfectMatches++;
							}
							else if(start.equals( Integer.toString(feature.location().bioStart())) ){
								startMatches++;
							}
							else if( end.equals( Integer.toString(feature.location().bioEnd()))){
								endMatches++;
							}
							else {
								notSameEndpointsMatches++;
							}
					}
				}
				firstIteration=false;
			}
			
			float freqPerfects = (float)perfectMatches/featuresToCompareTo.size() ;
			float freqStarts = (float)startMatches/featuresToCompareTo.size() ;
			float freqEnds = (float)endMatches/featuresToCompareTo.size() ;
			float freqNoEndpoints = (float)notSameEndpointsMatches/featuresToCompareTo.size() ;
			
			System.out.println("Fractions of annotated genes : Perfect, Start, End, NotSameEndpoints =");
			System.out.println(freqPerfects + " "+ freqStarts +" "+ freqEnds + " "+ freqNoEndpoints);
			
			 freqPerfects = (float)perfectMatches/NbGenesInFile ;
			 freqStarts = (float)startMatches/NbGenesInFile ;
			 freqEnds = (float)endMatches/NbGenesInFile ;
			 freqNoEndpoints = (float)notSameEndpointsMatches/NbGenesInFile ;
			 
			 System.out.println("Fractions of my predicted genes : Perfect, Start, End, NotSameEndpoints =");
			System.out.println(freqPerfects + " "+ freqStarts +" "+ freqEnds + " "+ freqNoEndpoints);
				
			
			fr.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		} 

	}

}
