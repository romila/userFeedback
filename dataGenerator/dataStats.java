import java.io.*;
import java.util.*;

public class dataStats {
	public List<List<String>> dataTuples;
	public int numSources;
	public int numObjects;
	public List<List<String>> objectValues;
	public int[] numberSourceVoted;
	
	public dataStats(String dataFile) throws IOException {
		this.dataTuples   = getTuplesFromFile (dataFile);
	
		this.numObjects   = this.dataTuples.size();
		this.numSources   = getNumberOfSources(this.dataTuples);
		this.numberSourceVoted = new int[this.numSources];
		getObjectUniqueValues();
	}
	
	/*
	 * function to read data file into tuples
	 */
	public List<List<String>> getTuplesFromFile(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		
		List<List<String>> fileTuples = new ArrayList<List<String>>();
        
    	String line = reader.readLine();
        while (line != null) {
        	String[] attributes = line.split("\t");
            List<String> singleTuple = new ArrayList<String>();
            for(int i=0; i<attributes.length; i++){
            	if(!attributes[i].equals(""))
            		singleTuple.add(attributes[i]);
            	else 
            		singleTuple.add(null);
            }
            line = reader.readLine();
            fileTuples.add(singleTuple);
        }        
        return fileTuples;
	}
	
	/*
	 * function to return number of sources from data tuples
	 */
	public int getNumberOfSources(List<List<String>> dataTuples) {
		List<Integer> numSourcesObjects = new ArrayList<Integer>();
		for (int i = 0; i < this.numObjects; i++)
			numSourcesObjects.add(dataTuples.get(i).size());
		return Collections.max(numSourcesObjects);
	}
	
	/*
	 * function to return unique object values for a tuple
	 */
	public void getObjectUniqueValues() {
		this.objectValues = new ArrayList<List<String>>();
		for (int i = 0; i < this.numObjects; i++) {
			ArrayList<String> objectUniqueValues = new ArrayList<String>();
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null) {
					objectUniqueValues.add(this.dataTuples.get(i).get(j));
					this.numberSourceVoted[j]++;
				}
			}
			objectUniqueValues = new ArrayList<String>(new HashSet<String>(objectUniqueValues));
			this.objectValues.add(objectUniqueValues);
		}
	}
	
	/*
	 * function to write results to file
	 */
	public  void writeToFile(double[] results, String fileName) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
		
		for ( int i = 0; i < results.length; i++)
	        writer.write(results[i] + "\n");
		
		writer.close();
	}

	/*
	 * function to get degree of each object
	 */
	public double[] getDegreeOfObjects() {
		double[] count = new double[this.dataTuples.size()];
		for (int i = 0; i < this.dataTuples.size(); i++) {
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null)
					count[i]++;
			}
		}
		return count;
	}
	
	/* 
	 * function to get agreement/disagreement for each object,
	 * levelOfConflict .. ratio of conflict in top contesting values
	 */
	public double[] getConflictRatios() {
		double[] conflictRatio = new double[this.dataTuples.size()];
		for (int i = 0; i < this.dataTuples.size(); i++) {
			List<String> objectValues = this.objectValues.get(i);
			objectValues.removeAll(Collections.singleton(null));
			
			String[][] valuesCount = new String[objectValues.size()][2]; // value in [0], count in [1]
			for (int j = 0; j < valuesCount.length; j++) {
				valuesCount[j][0] = objectValues.get(j);
				valuesCount[j][1] = String.valueOf(0);
			}
			
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null) {
					for (int k = 0; k < valuesCount.length; k++) {
						if (valuesCount[k][0].equals(this.dataTuples.get(i).get(j))) {
							valuesCount[k][1] = String.valueOf(Integer.parseInt(valuesCount[k][1]) + 1);
							break;
						}
					}
				}
			}
			
			if (objectValues.size() > 0) {
				int firstMax = 0; int index = 0;
				int secondMax = 0;
				
				for (int j = 0; j < valuesCount.length; j++) {
					if (Integer.parseInt(valuesCount[j][1]) > firstMax) {
						firstMax = Integer.parseInt(valuesCount[j][1]);
						index = j;
					}
				}
				
				valuesCount[index][1] = String.valueOf(0); index = 0;
				for (int j = 0; j < valuesCount.length; j++) {
					if (Integer.parseInt(valuesCount[j][1]) > secondMax) {
						secondMax = Integer.parseInt(valuesCount[j][1]);
						index = j;
					}
				}
				
				if (secondMax != 0)
					conflictRatio[i] = (double) firstMax/secondMax;
			}
		}
		
		return conflictRatio;
	}
	
	/*
	 * function to get entropy for each object
	 */
	public double[] getEntropyOfObjects() {
		double[] entropy = new double[this.dataTuples.size()];
		
		for (int i = 0; i < this.dataTuples.size(); i++) {
			int count = 0;
			for (int j = 0; j < this.dataTuples.get(i).size(); j++)
				if (this.dataTuples.get(i).get(j) != null)
					count++;
			
			List<String> objectValues = this.objectValues.get(i);
			objectValues.removeAll(Collections.singleton(null));
			
			if (objectValues.size() > 0) {
				double[] fraction = new double[objectValues.size()];
				for (int j = 0; j < objectValues.size(); j++) {
					fraction[j] = (double) Collections.frequency(this.dataTuples.get(i), objectValues.get(j)) / count;
					if (fraction[j] > 0)
						entropy[i] += -1 * fraction[j] * Math.log(fraction[j]); 
				}
			}
		}
		
		return entropy;
	}
	
	/*
	 * function to get gini value for each object
	 */
	public double[] getGiniValueOfObjects() {
		double[] giniValue = new double[this.dataTuples.size()];
		
		for (int i = 0; i < this.dataTuples.size(); i++) {
			int count = 0;
			for (int j = 0; j < this.dataTuples.get(i).size(); j++)
				if (this.dataTuples.get(i).get(j) != null)
					count++;
			
			List<String> objectValues = this.objectValues.get(i);
			objectValues.removeAll(Collections.singleton(null));
			
			if (objectValues.size() > 0) {
				double[] fraction = new double[objectValues.size()];
				for (int j = 0; j < objectValues.size(); j++) {
					fraction[j] = (double) Collections.frequency(this.dataTuples.get(i), objectValues.get(j)) / count;
					if (fraction[j] > 0)
						giniValue[i] += fraction[j] * fraction[j]; 
				}
			}
			giniValue[i] = 1 - giniValue[i];
		}
		return giniValue;
	}
	
	/*
	 * function to get % of sources that vote for an object
	 */
	public double[] getFractionOfSourcesVotingForObject() {
		double[] fractionVoted = new double[this.dataTuples.size()];
		for (int i = 0; i < this.dataTuples.size(); i++) {
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null)
					fractionVoted[i]++;
			}
			fractionVoted[i] /= this.numSources;
		}
		return fractionVoted;
	}
	
	/*
	 *  args[0]: filename that needs to be analyzed
	 */
	public static void main(String[] args) throws IOException {
		dataStats d = new dataStats(args[0]);
		
		double[] degreeOfObjects = d.getDegreeOfObjects();
		double[] conflictOfObjects = d.getConflictRatios();		
		double[] entropyOfObjects = d.getEntropyOfObjects();
		double[] giniValueOfObjects = d.getGiniValueOfObjects();
		double[] fractionOfSourcesVotingForObjects = d.getFractionOfSourcesVotingForObject();
		
		d.writeToFile(degreeOfObjects, "degree.txt");
		d.writeToFile(conflictOfObjects, "conflictRatio.txt");
		d.writeToFile(entropyOfObjects, "entropy.txt");
		d.writeToFile(giniValueOfObjects, "giniValue.txt");
		d.writeToFile(fractionOfSourcesVotingForObjects, "fractionVoting.txt");
		
		System.out.println("End of program");
	}
}
