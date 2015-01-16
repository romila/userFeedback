import java.io.*;
import java.util.*;

public class generator {
	public String degreeDistributionFile;
	public String giniDistributionFile;
	public String fractionVotingFile;
	public String conflictRatioFile;
	
	public int numObjects; // #objects in synthetic data
	public int numSources; // #sources in synthetic data
	public String data[][];
	public String truth[][];
		
	public generator(String[] args) {
		this.numObjects = Integer.parseInt(args[0]);
		this.numSources = Integer.parseInt(args[1]);
		
		this.degreeDistributionFile = args[2];
		this.giniDistributionFile = args[3]; 
		this.fractionVotingFile = args[4];
		this.conflictRatioFile = args[5];
		
		this.data = new String[this.numObjects][this.numSources];
		this.truth = new String[this.numObjects][2];
	}
	
	/*
	 * function to read a file
	 */
	public List<String> getValuesFromFile(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
			
		List<String> fileTuples = new ArrayList<String>();
        
    	String line = reader.readLine();
        while (line != null) {
        	fileTuples.add(line);
        	line = reader.readLine();
        }
        reader.close();
        return fileTuples;
	}
	
	/*
	 * function to write results to a file
	 */
	public  void writeToFile(String[][] results, String fileName) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			
		for ( int i = 0; i < results.length; i++) {
			for(int j = 0; j < results[i].length; j++ ) {
				if (results[i][j] != null)
					writer.write(results[i][j]);
				writer.write("\t");
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	/*
	 * function to return unique object values for a tuple
	 */
	public List<String> getObjectUniqueValues(String[] values) {
		ArrayList<String> uniqueValues = new ArrayList<String>();
		for (int i = 0; i < values.length; i++) {
			if (values[i] != null)
				uniqueValues.add(values[i]);
		}
		
		return (new ArrayList<String>(new HashSet<String>(uniqueValues)));
	}
	
	/* 
	 * function to write to data and truth files
	 */
	public void generateDataAndTruthFiles() {
		try {
			writeToFile(this.data, "generatedData_" + this.numObjects + ".txt");
		} catch (IOException e) {
			System.out.println("Writing error: data file");
		}
		
		try {
			writeToFile(this.truth, "generatedTruth_" + this.numObjects + ".txt");
		} catch (IOException e) {
			System.out.println("Writing error: data file");
		}
	}
	
	/*
	 * function to generate completely random data and truth
	 */
	public void randomData() {
		String[] values = {"0", "1"};
		int randomValue = 0;
		
		Random random = new Random(System.currentTimeMillis());
		for (int i = 0; i < this.numObjects; i++) {
			for (int j = 0; j < this.numSources; j++) {
				randomValue = random.nextInt(3);
				if (randomValue != 2)
					this.data[i][j] = values[randomValue];
			}
			
			List<String> uniqueValues = this.getObjectUniqueValues(this.data[i]);
			if (uniqueValues.size() == 0) { // no source reported a value
				Random r = new Random(System.currentTimeMillis());
				this.data[i][r.nextInt(this.numSources)] = values[0]; // have a random source report a value
				
				this.truth[i][0] = String.valueOf(i);
				this.truth[i][1] = String.valueOf(values[0]);
			}
			
			else {
				this.truth[i][0] = String.valueOf(i);
				this.truth[i][1] = String.valueOf(uniqueValues.get(random.nextInt(uniqueValues.size())));
			}
		}
	}
	
	/*
	 * function to generate data with similar distribution for #votes
	 */
	public void degreeData() throws IOException {
		List<String> oldDegrees = this.getValuesFromFile(this.degreeDistributionFile);
		
		double realMean = 0;
		double realStddev = 0;
		
		for (int i = 0; i < oldDegrees.size(); i++) 
			realMean += Double.parseDouble(oldDegrees.get(i));
		realMean /= oldDegrees.size();
		
		for (int i = 0; i < oldDegrees.size(); i++) 
			realStddev += Math.pow((Double.parseDouble(oldDegrees.get(i)) - realMean), 2);
		realStddev = Math.sqrt(realStddev/(double) oldDegrees.size());
		
		// adjusting mean and stddev for degrees
		double newMean = ((double) this.numObjects / oldDegrees.size()) * realMean;
		double newStddev = ((double) this.numObjects / oldDegrees.size()) * realStddev;
		
		String[] values = {"0", "1"};
		double newDegree;
		Random r = new Random(System.currentTimeMillis()); // for normal distribution generator
		for (int i = 0; i < this.numObjects; i++) {
			newDegree = Math.min(Math.ceil(r.nextGaussian() * newStddev + newMean), this.numSources);
			for (int j = 0; j < newDegree; j++) {
				Random r0 = new Random(System.currentTimeMillis()); // for location of vote
				Random r1 = new Random(System.currentTimeMillis()); // for value of vote
				int vote = r1.nextInt(3);
				if (vote != 2)
					this.data[i][r0.nextInt(this.numSources)] = values[vote];
			}
			
			List<String> uniqueValues = this.getObjectUniqueValues(this.data[i]);
			if (uniqueValues.size() == 0) { // no source reported a value
				Random r2 = new Random(System.currentTimeMillis());
				this.data[i][r2.nextInt(this.numSources)] = values[0]; // have a random source report a value
				
				this.truth[i][0] = String.valueOf(i);
				this.truth[i][1] = String.valueOf(values[0]);
			}
			
			else {
				Random r3 = new Random(System.currentTimeMillis());
				this.truth[i][0] = String.valueOf(i);
				this.truth[i][1] = String.valueOf(uniqueValues.get(r3.nextInt(uniqueValues.size())));
			}
		}
	}
		
	/* args[0]: n (#objects)
	 * args[1]: #sources
	 * args[2]: degreeFile,
	 * args[3]: giniFile
	 * args[4]: fractionVotingFile
	 * args[5]: conflictRatioFile
	 */
	public static void main(String[] args) throws IOException {
		generator g = new generator(args);
		g.randomData();
		g.generateDataAndTruthFiles();

		g.degreeData();
		g.generateDataAndTruthFiles();

		
	System.out.println("End \n");
	}
}
