import java.io.*;
import java.util.*;

public class generator {
	public String degreeDistributionFile;
	public String entropyDistributionFile;
	public int numObjects;
	public int numSources;
	public String data[][];
	public String truth[][];
	public Random random; 
		
	public generator (int nObjects, int nSources) {
		this.numObjects = nObjects;
		this.numSources = nSources;
		this.random = new Random(System.currentTimeMillis());
		this.data = new String[this.numObjects][this.numSources];
		this.truth = new String[this.numObjects][2];
	}
	
	// read file
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
	
	// write results to file
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
	
	// generate nObjects with normal distribution
	public double[] drawFromNormalDistribution(List<String> distributionStrings, int nObjects) {
		double mean = 0;
		double stddev = 0;
		Random random = new Random(System.currentTimeMillis());
		
		for (int i = 0; i < distributionStrings.size(); i++) 
			mean += Double.parseDouble(distributionStrings.get(i));
		
		mean /= distributionStrings.size();
		
		for (int i = 0; i < distributionStrings.size(); i++) 
			stddev += Math.pow((Double.parseDouble(distributionStrings.get(i)) - mean), 2);
		
		stddev = Math.sqrt(stddev/(double) distributionStrings.size());
		
		// adjusting mean and stddev for degrees
//		if (mean > 1 && stddev > 1) {
//			double oldMean = mean;
//			mean *= this.numSources/(double)38; // 38 sources in real data
//			stddev *= mean/oldMean;
//		}
		
		double[] values = new double[nObjects];
		int count = 0;

		// adjusting mean and stddev for degrees
		if (mean > 1 && stddev > 1) {
			stddev *= (this.numSources/2)/mean;
			mean = Math.ceil((double) this.numSources/2);
		
			double tempValue;
			boolean withinRange = false;
			while (count < nObjects) {
				withinRange = false;
				while (!withinRange) {
					tempValue = Math.abs(random.nextGaussian() * stddev + mean); // degree normal distribution
					if (tempValue <= this.numSources && tempValue > 1) { //at least two sources and at most n sources vote
						values[count] = tempValue;
						count++;
						withinRange = true;
					}
				}
			}
		}
		
		else {
			count = 0;
			while (count < nObjects) {
				values[count] = Math.abs(random.nextGaussian() * stddev + mean);
				count++;
			}
		}
		
		return values;
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
	
	/* args[0]: n (#objects)
	 * args[1]: #sources
	 * args[2]: degreeFile,
	 * args[3]: entropyFile
	 */
	public static void main(String[] args) {
		int nObjects = Integer.parseInt(args[0]);
		int nSources = Integer.parseInt(args[1]);
		
		generator g = new generator(nObjects, nSources);
		List<String> uniqueValues = null;
		
		if (args.length == 2) { // random data
			String[] values = {"0", "1"};
			
			for (int i = 0; i < g.numObjects; i++) {
				for (int j = 0; j < g.numSources; j++) {
					int randomValue = g.random.nextInt(3);
					if (randomValue != 2)
						g.data[i][j] = values[randomValue];
				}
				uniqueValues = g.getObjectUniqueValues(g.data[i]);
				if (uniqueValues.size() == 0)
					g.data[i][0] = values[0];
				uniqueValues = g.getObjectUniqueValues(g.data[i]);
				g.truth[i][0] = String.valueOf(i);
				g.truth[i][1] = String.valueOf(uniqueValues.get(g.random.nextInt(uniqueValues.size()))); 
			}
		}
		
		if (args.length == 4) { // data follows degree and entropy distribution
			g.degreeDistributionFile = args[2];
			g.entropyDistributionFile = args[3];
			
			List<String> degreeStrings = new ArrayList<String>();
			List<String> entropyStrings = new ArrayList<String>();
			
			try {
				degreeStrings = g.getValuesFromFile(args[2]);
			} catch (IOException e) {
				System.out.println("Reading error: degree file");
			}
			
			try {
				entropyStrings = g.getValuesFromFile(args[3]);
				entropyStrings.removeAll(Collections.singleton("0.0"));
			} catch (IOException e) {
				System.out.println("Reading error: entropy file");
			} 
			
			double[] degreeOfObjects = g.drawFromNormalDistribution(degreeStrings, g.numObjects);
			for (int i = 0; i < degreeOfObjects.length; i++) 
				degreeOfObjects[i] = Math.ceil(degreeOfObjects[i]);
			
			// barring single-valued objects, entropy follows a normal distribution
			double[] entropyOfObjects = g.drawFromNormalDistribution(entropyStrings, g.numObjects);
			
			int countDegree = 0;
			int countEntropy = 0;
			double differentVotes = 0;
			String[] values = {"0", "1"};
			
			for (int i = 0; i < g.numObjects; i++) {
				boolean[] votedBy = new boolean[g.numSources];
				countDegree = 0;
				countEntropy = 0;
				System.out.println(i + ":\t" + degreeOfObjects[i] + "\t" + entropyOfObjects[i]);
				if (entropyOfObjects[i] <= 0.5) { // real roots
					differentVotes = Math.ceil((1 - Math.sqrt(1 - 2 * entropyOfObjects[i]))/2 * degreeOfObjects[i]);
					if(Double.isNaN(differentVotes))
						differentVotes = Math.ceil((1 + Math.sqrt(1 - 2 * entropyOfObjects[i]))/2 * degreeOfObjects[i]);
				}
				else // imaginary roots
					differentVotes = Math.ceil(degreeOfObjects[i]/2);
				
				// take degree distribution into account 
				while (countDegree < degreeOfObjects[i]) {
					int randomSource = g.random.nextInt(g.numSources);
					if (!votedBy[randomSource]) {
						votedBy[randomSource] = true;
						countDegree++;
						// take entropy distribution into account
						if (countEntropy < differentVotes) {
							g.data[i][randomSource] = String.valueOf(1);
							countEntropy++;
						}
						else 
							g.data[i][randomSource] = String.valueOf(0);
				 	}
				}
				uniqueValues = g.getObjectUniqueValues(g.data[i]);
				if (uniqueValues.size() == 0)
					g.data[i][0] = values[0];
				
				g.truth[i][0] = String.valueOf(i);
				g.truth[i][1] = String.valueOf(uniqueValues.get(g.random.nextInt(uniqueValues.size())));
			}
		}
		
		g.generateDataAndTruthFiles();

		
	System.out.println("End \n");
	}
}
