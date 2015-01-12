import java.io.*;
import java.util.*;

public class dataStats {
	
	// read data file into tuples
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
	
	// write results to file
	public  void writeToFile(List<String> results, String fileName) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
		
		for ( int i = 0; i < results.size(); i++)
	        writer.write(results.get(i)+ "\n");
		
		writer.close();
	}

	// get degree of each object
	public int getDegreeOfObject(List<String> tupleData) {
		int count = 0;
		for (int i = 0; i < tupleData.size(); i++) 
			if (tupleData.get(i) != null)
				count++;
		
		return count;
	}
	
	/* get agreement/disagreement for each object,
	 * levelOfConflict .. ratio of conflict in top contesting values
	 */
	public double getConflictRatio(List<String> tupleData) {
		if (tupleData.size() == 0)
			return 0;
		
		ArrayList<String> objectUniqueValues = new ArrayList<String>(new HashSet<String>(tupleData));
		objectUniqueValues.removeAll(Collections.singleton(null));
		
		String[][] valuesCount = new String[objectUniqueValues.size()][2]; // value in [0], count in [1]
		for (int i = 0; i < valuesCount.length; i++) {
			valuesCount[i][0] = objectUniqueValues.get(i);
			valuesCount[i][1] = String.valueOf(0);
		}
		
		for (int i = 0; i < tupleData.size(); i++) {
			if (tupleData.get(i) != null) {
				for (int j = 0; j < valuesCount.length; j++) {
					if (valuesCount[j][0].equals(tupleData.get(i))) {
						valuesCount[j][1] = String.valueOf(Integer.parseInt(valuesCount[j][1]) + 1);
						break;
					}
				}
			}
		}
		
		int firstMax = 0; int index = 0;
		int secondMax = 0;
		
		for (int i = 0; i < valuesCount.length; i++) {
			if (Integer.parseInt(valuesCount[i][1]) > firstMax) {
				firstMax = Integer.parseInt(valuesCount[i][1]);
				index = i;
			}
		}
		
		valuesCount[index][1] = String.valueOf(0); index = 0;
		for (int i = 0; i < valuesCount.length; i++) {
			if (Integer.parseInt(valuesCount[i][1]) > secondMax) {
				secondMax = Integer.parseInt(valuesCount[i][1]);
				index = i;
			}
		}
		
		double conflictRatio = 0;
		if (secondMax != 0)
			conflictRatio = (double) firstMax/secondMax;
		System.out.println(conflictRatio);
		
		return conflictRatio;
		
		
	}
	
	// get entropy of each object
	public double getEntropyOfObject(List<String> tupleData) {
		double entropy = 0;
		int count = 0;
		for (int i = 0; i < tupleData.size(); i++) 
			if (tupleData.get(i) != null)
				count++;
		
		ArrayList<String> objectUniqueValues = new ArrayList<String>(new HashSet<String>(tupleData));
		objectUniqueValues.removeAll(Collections.singleton(null));
		
		if (objectUniqueValues.size() > 0) {
			double[] fractionUnique = new double[objectUniqueValues.size()];
			for (int i = 0; i < objectUniqueValues.size(); i++) {
//				System.out.println(Collections.frequency(tupleData, objectUniqueValues.get(i)));
				fractionUnique[i] = (double) Collections.frequency(tupleData, objectUniqueValues.get(i))/count;
				if (fractionUnique[i] > 0)
					entropy += -1 * fractionUnique[i] * Math.log(fractionUnique[i]);
			}
		}
		
		return entropy;
	}
	
	// get entropy of each object
	public double getGiniValueOfObjects(List<String> tupleData) {
		double giniValue = 0;
		int count = 0;
		for (int i = 0; i < tupleData.size(); i++) 
			if (tupleData.get(i) != null)
				count++;
		
		ArrayList<String> objectUniqueValues = new ArrayList<String>(new HashSet<String>(tupleData));
		objectUniqueValues.removeAll(Collections.singleton(null));
		
		if (objectUniqueValues.size() > 0) {
			double[] fractionUnique = new double[objectUniqueValues.size()];
			for (int i = 0; i < objectUniqueValues.size(); i++) {
				fractionUnique[i] = (double) Collections.frequency(tupleData, objectUniqueValues.get(i))/count;
				if (fractionUnique[i] > 0)
					giniValue += fractionUnique[i] * fractionUnique[i];
			}
		}
		
		return (1 - giniValue);
	}
	
	// args[0]: filename
	public static void main(String[] args) {
		dataStats d = new dataStats();
		List<List<String>>  dataTuples = new ArrayList<List<String>>();
		try {
			dataTuples = d.getTuplesFromFile (args[0]);
		}
		catch (Exception e) {
			System.out.println("File reading error");
		}
		
		List<String> degreeOfObjects = new ArrayList<String>();
		for (int i = 0; i < dataTuples.size(); i++)
			degreeOfObjects.add(String.valueOf(d.getDegreeOfObject(dataTuples.get(i))));

		try {
			d.writeToFile(degreeOfObjects, "degreeOfObjects.txt");
		} catch (Exception e) {
			System.out.println("File writing error");
		}
		
//		List<Double> conflictOfObjects = new ArrayList<Double>();
//		for (int i = 0; i < dataTuples.size(); i++)
//			conflictOfObjects.add(d.getConflictRatio(dataTuples.get(i)));
//
//		try {
//			d.writeToFile(conflictOfObjects, "conflictOfObjects.txt");
//		} catch (Exception e) {
//			System.out.println("File writing error");
//		}

//		List<Double> entropyOfObjects = new ArrayList<Double>();
//		for (int i = 0; i < dataTuples.size(); i++)
//			entropyOfObjects.add(d.getEntropyOfObject(dataTuples.get(i)));
//
//		try {
//			d.writeToFile(entropyOfObjects, "entropyOfObjects.txt");
//		} catch (Exception e) {
//			System.out.println("File writing error");
//		}
		
		List<String> giniValueOfObjects = new ArrayList<String>();
		for (int i = 0; i < dataTuples.size(); i++)
			giniValueOfObjects.add(String.valueOf(d.getGiniValueOfObjects(dataTuples.get(i))));

		try {
			d.writeToFile(giniValueOfObjects, "giniValueOfObjects.txt");
		} catch (Exception e) {
			System.out.println("File writing error");
		}
	}
}
