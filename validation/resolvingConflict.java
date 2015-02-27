import java.io.*;
import java.util.*;
import java.nio.file.Paths;

/*
 * ASSUMPTIONS:
 * dataFile:
 *  - one object is voted by more than one source
 *  - each object has at least two different votes
 *  - exactly one of all the votes is true, rest are false
 *  
 * truthFile:
 *  - true values are known for a sample of the dataFile
 *  - file might not have objects in the same order as in the dataFile
 */

public class resolvingConflict {
	public List<List<String>> dataTuples;
	public List<List<String>> truthTuples;
	public int numObjects;
	public int numSources;
	public List<List<String>> objectValues; // List of (List of different values for each object)
	public int[] numberSourceVoted; // array containing the number of objects a source voted for
	public accuPR basePredictor;
	
	public resolvingConflict(String dataFile, String truthFile) throws IOException {
		this.dataTuples   = getTuplesFromFile (dataFile);
		this.truthTuples  = getTuplesFromFile (truthFile);
		
		this.numObjects   = this.dataTuples.size();
		this.numSources   = getNumberOfSources(this.dataTuples);
		
		this.numberSourceVoted = new int[this.numSources];
		getObjectUniqueValues();
		this.basePredictor = new accuPR (this.dataTuples, null);
	}
	
	/*
	 * function to read tuples from a file
	 */
	public List<List<String>> getTuplesFromFile(String fileName) throws IOException {
		try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
			List<List<String>> fileTuples = new ArrayList<List<String>>();
		    String line = reader.readLine();
	        while (line != null) {
	        	String[] attributes = line.split("\t");
	            List<String> singleTuple = new ArrayList<String>();
	            for (int i = 0; i < attributes.length; i++) {
	            	if (!attributes[i].equals("")) {
	            		singleTuple.add(attributes[i]);
	            	}
	            	else 
	            		singleTuple.add(null);
	            }
	            line = reader.readLine();
	            fileTuples.add(singleTuple);
	        }       
	        reader.close();
			return fileTuples;
		}
	}
	
	/*
	 * function to write a 2D array of values to a file
	 */
	public void write_2D_ToFile(double[][] objectsGraph, String fileName) throws IOException {
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
			for (int i = 0; i < objectsGraph.length; i++) {
				for (int j = 0; j < objectsGraph[i].length; j++)
					writer.write(objectsGraph[i][j] + "\t");
				writer.write("\n");
			}
			writer.close();
		}
	}
	
	/*
	 * function to write a 1D array of values to a file
	 */
	public void write_1D_ToFile(double[] results, String fileName) throws IOException {
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
			for (int i = 0; i < results.length; i++) 
					writer.write(results[i] + "\n");
			writer.close();
		}
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
	 * function to return index having maximum in a set of values
	 */
	public int getIndexOfMaximumValue(double[] values) {
		double max = Integer.MIN_VALUE; int index = 0; 
		for (int i = 0; i < values.length; i++) {
			if (values[i] > max) {
					max = values[i];
					index = i;
			}
		}
		return index;
	}
	
	/*
	 * function to sort an array in decreasing order
	 * returns locations of sorted items in actual array
	 * e.g. [4, 5, 3, 6] returns [3, 1, 0, 2]
	 */
	public int[] sortAndRank(double[] a) {
		int j;
		double tempVal;
		int tempRank;
		
		int[] rank = new int[a.length];

		for (int i = 1; i < a.length; i++) {
			j = i;
			rank[i] = i;

			while (j > 0 && a[j-1] < a[j]) {
				tempVal = a[j-1];
				a[j-1] = a[j];
				a[j] = tempVal;

				tempRank = rank[j-1];
				rank[j-1] = rank[j];
				rank[j] = tempRank;

				j--;
			}
		}
		return rank;
	}
	
	/* 
	 * function to compute effectiveness assuming groundTruth is known for the complete database
	 * computing distance from truth 
	 * distanceFromTruth d = average over all objects (1 - probability of true value)
	 * *myAngle: d* for one object = sum over all distances from absolute probability distribution
	 * *i.e., if truth = {1, 0, 0} and p = {0.7, 0.2, 0.1}
	 * d* = |1 - 0.7| + |0 - 0.2| + |0 - 0.1| = 0.6, while
	 * d =  |1 - 0.7| = 0.3  
	 * conceptually, d* covers the entire distance but mathematically, while comparing two objects,
	 * d and d* differ only by a factor of 2.
	 */
	public double computeDistanceFromTruth(accuPR t) {
		double distance = 0;
		int currentObject;
		int locationOfTrueValue = 0;
		List<String> currentObjectValues = new ArrayList<String>();
		
		for (int i = 0; i < this.truthTuples.size(); i++) {
			currentObject = Integer.parseInt(this.truthTuples.get(i).get(0));
			currentObjectValues = this.objectValues.get(currentObject);
			if (this.truthTuples.get(i).size() > 1) {
				locationOfTrueValue = currentObjectValues.indexOf(this.truthTuples.get(i).get(1));
				if (locationOfTrueValue == -1) // if none of the object values is true
					distance += 1;
				else
					distance += Math.abs(1 - t.getValueProbability()[currentObject][locationOfTrueValue]);
			}
		}
		distance /= this.truthTuples.size();
		return distance;
	}
	
	/*
	 * function to compute overall database uncertainty 
	 * computed as Shannon entropy, U_d = sum over all objects {E}
	 * where E = sum over all values {-p log(p)}
	 */
	public double computeDatabaseUncertainty(double[][] probabilityValues) {
		double dbUncertainty = 0; 
		for (int i = 0; i < probabilityValues.length; i++) {
			for (int j = 0; j < probabilityValues[i].length; j++) {
				if (probabilityValues[i][j] > 0)
					dbUncertainty += -1 * probabilityValues[i][j] * Math.log(probabilityValues[i][j]);
			}
		}
		return dbUncertainty;
	}
		
	/*
	 * function to select object that results in maximum database utility gain (using ground truth)
	 * MUG: Maximum Utility Gain
	 */
	public int selectObjectWithMUG(List<Integer> tuplesNotValidated, List<List<String>> indices) {
		
		double[] listDistance 	 = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++) {
			List<String> truthTuple = new ArrayList<String>(this.truthTuples.get(tuplesNotValidated.get(i)));
			List<List<String>> tempIndices = new ArrayList<List<String>>(indices);
			tempIndices.add(truthTuple);
			accuPR	tempValidatedPredictor	= new accuPR(this.dataTuples, tempIndices);
			listDistance[i] = 1 - computeDistanceFromTruth(tempValidatedPredictor); 
			// distance from truth: 7% < 18%, Gain: 1-7% > 1-18%
		}
		
		// get object that is closest to truth now
		int indexMax = getIndexOfMaximumValue(listDistance);
		return indexMax;
	}

	/*
	 * function to select object that results in maximum expected database utility gain
	 * MEUG: Maximum Expected Utility Gain
	 */
	public int selectObjectWithMEUG(List<Integer> tuplesNotValidated, List<List<String>> indices) {
		
		double baseDatabaseUncertainty 	= computeDatabaseUncertainty (this.basePredictor.getValueProbability());
		double[ ] expectedUtility 	 = new double[tuplesNotValidated.size()];
		
		List<String> indexTuple = new ArrayList<String>();
		for (int i = 0; i < tuplesNotValidated.size(); i++) {
			double tupleDBUncertainty = 0;
			double denominator = 0;
			for (int j = 0; j < this.objectValues.get(tuplesNotValidated.get(i)).size(); j++) {
				List<List<String>> tempIndices = new ArrayList<List<String>>(indices);
				
				indexTuple.clear();
				indexTuple.add(Integer.toString(tuplesNotValidated.get(i))); 
				indexTuple.add(this.objectValues.get(tuplesNotValidated.get(i)).get(j));
				
				tempIndices.add(indexTuple);
				
				accuPR	tempValidatedPredictor	= new accuPR(this.dataTuples, tempIndices);
				tupleDBUncertainty += this.basePredictor.getValueProbability()[tuplesNotValidated.get(i)][j] * 
										computeDatabaseUncertainty(tempValidatedPredictor.getValueProbability());
				denominator += this.basePredictor.getValueProbability()[tuplesNotValidated.get(i)][j];
			}
			tupleDBUncertainty /= denominator;
			
			expectedUtility[i] = - tupleDBUncertainty + baseDatabaseUncertainty; 
			// (tupUncer - baseUncer) should be -ve. Meaning uncertainty has been decreased. Thus, look for object
			// that has minimum(tupUncer - baseUncer) or maximum(-tupUncer + baseUncer)
		}
		
		int indexMax = getIndexOfMaximumValue(expectedUtility); // object that results in maximum utility
		return indexMax;
	}
	
	/*
	 * function to compute validation by selecting object that results in maximum database utility gain
	 * using ground truth
	 */
	public double[] computeForMU() throws IOException {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		double[] distancesComputed = new double[this.numObjects + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		int count = 0;
		
		accuPR tempPredictor;
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		
		int countValidated = 0;
		int countValidatedIterator = 1;
		
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < 50 * countValidatedIterator &&
					!tuplesNotValidated.isEmpty()) {
		
				int index = selectObjectWithMUG(tuplesNotValidated, indices);
				indices.add(this.truthTuples.get(tuplesNotValidated.get(index)));
				tuplesNotValidated.remove(index);
				countValidated++;	
			}
				
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
		
			if (indices.size() > 100 && indices.size() % ((indices.size()/100) * 100) >= 0)
				System.out.println(indices.size() + " : " + countValidated);
			countValidatedIterator++;
			
		}
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object that results in maximum expected database utility gain
	 */
	public double[]  computeForMEU() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		double[] distancesComputed = new double[this.numObjects + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		int count = 0;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		accuPR tempPredictor;
		
		int countValidated = 0;
		int countValidatedIterator = 1;
		
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < 50 * countValidatedIterator &&
					!tuplesNotValidated.isEmpty()) {
				int index = selectObjectWithMEUG(tuplesNotValidated, indices);
				indices.add(this.truthTuples.get(tuplesNotValidated.get(index)));
				tuplesNotValidated.remove(index);
				countValidated++;	
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			
			if (indices.size() > 100 && indices.size() % ((indices.size()/100) * 100) >= 0)
				System.out.println(indices.size() + " : " + countValidated);
			countValidatedIterator++;
		}
				
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting objects randomly
	 */
	public double[] computeForRandomSelection(int numberOfRuns) {
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // list of <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		Random random = new Random(System.currentTimeMillis());
		accuPR tempPredictor;
	    
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		
		for (int i = 0; i < this.truthTuples.size(); i++) {
			if (this.truthTuples.get(i).size() > 1)
				tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		}
		
		for (int i = 0; i < numberOfRuns; i++) {
			count = 1;
			indices.clear();
			
			int countValidated = 0;
			int countValidatedIterator = 1;
			while (!tuplesNotValidated.isEmpty()) {
				while (indices.size() < 50 * countValidatedIterator  && 
						!tuplesNotValidated.isEmpty()) {
					int index = random.nextInt(tuplesNotValidated.size());
				    indices.add(this.truthTuples.get(tuplesNotValidated.get(index)));
					tuplesNotValidated.remove(index);
					countValidated++;				
				}
				
				tempPredictor = new accuPR(this.dataTuples, indices); 
				distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
				
				if (indices.size() > 100 && indices.size() % ((indices.size()/100) * 100) >= 0)
					System.out.println(indices.size() + " : " + countValidated);
				countValidatedIterator++;
			}
					
			// again add all objects to unvalidated list for next run
			for (int j = 0; j < this.truthTuples.size(); j++) {
				if (this.truthTuples.get(j).size() > 1)
					tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(j).get(0)));
			}
		}
				
		for (int i = 1; i < distancesComputed.length; i++)
			distancesComputed[i] /= numberOfRuns;
		
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object based on their entropies
	 * aka MajorityVotingOrdering, or MVO
	 */
	public double[] computeForMVO() {
		double fraction = 0;
		double countNotNull = 0;
		int currentObject;
		double[] objectEntropies = new double[this.truthTuples.size()];
		List<String> currentValues = new ArrayList<String>();
		
		// get entropies of objects in the truthfile
		for (int i = 0; i < this.truthTuples.size(); i++) {
			currentObject = Integer.parseInt(this.truthTuples.get(i).get(0));
			List<String> objectTupleList = this.dataTuples.get(currentObject);
			countNotNull = objectTupleList.size() - Collections.frequency(objectTupleList, null);
			
			currentValues = this.objectValues.get(currentObject);
			
			if (currentValues.size() > 1) {
				for (int j = 0; j < currentValues.size(); j++) {
					fraction = (double) Collections.frequency(objectTupleList, currentValues.get(j)) / countNotNull;
					objectEntropies[i]  += -1 * fraction * Math.log(fraction);
				}
			}
		}
		
		// sort objects in decreasing order of their entropies
		int[] entropyRanks = sortAndRank(objectEntropies);
		
		accuPR tempPredictor;
		
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		int countValidated = 0;
		int countValidatedIterator = 1;
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < 50 * countValidatedIterator &&
					!tuplesNotValidated.isEmpty()) {
				outerloop:
				for (int i = countValidated; i < entropyRanks.length; i++) {
					if (tuplesNotValidated.contains(entropyRanks[i])) {
						indices.add(this.truthTuples.get(truthIndices.indexOf(entropyRanks[i])));
						tuplesNotValidated.remove(tuplesNotValidated.indexOf(entropyRanks[i]));
						countValidated++;
						break outerloop;
					}
				}
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			
			if (indices.size() > 100 && indices.size() % ((indices.size()/100) * 100) >= 0)
				System.out.println(indices.size() + " : " + countValidated);
			countValidatedIterator++;
		}
		return distancesComputed;
	}
	
	/*
	 * for 2-valued case, returns the cardinality of set of sources that vote for
	 * value i in A and j in B, {i,j}
	 * {{dA,dB}, {dA,nB}, {nA,dB}, {nA,nB}}
	 */
	public int[] numSources_set(List<String> aValues, List<String> bValues,
								String domA, String nonDomA,
								String domB, String nonDomB) {
		int[] num_set = new int[4]; 
		
		for (int i = 0; i < this.numSources; i++) {
			if (aValues.size() > i && bValues.size() >i &&
					aValues.get(i) != null && bValues.get(i) != null) {
				if (aValues.get(i).equals(domA) && bValues.get(i).equals(domB))
					num_set[0]++;
				else if (aValues.get(i).equals(domA) && bValues.get(i).equals(nonDomB))
					num_set[1]++;
				else if (aValues.get(i).equals(nonDomA) && bValues.get(i).equals(domB))
					num_set[2]++;
				else if (aValues.get(i).equals(nonDomA) && bValues.get(i).equals(nonDomB))
					num_set[3]++;
			}
		}
		return num_set;
	}
	
	/*
	 * function to compute impact on b due to validation of a
	 * 2 - value case
	 */
	public double getDelta_b_dueTo_a(int a, int b) {
		List<String> aValues = this.dataTuples.get(a);
		List<String> bValues = this.dataTuples.get(b);
		List<String> a_uniqueValues = this.objectValues.get(a);
		List<String> b_uniqueValues = this.objectValues.get(b);
		
		if (aValues.size() == 0 || bValues.size() == 0)
			return 0;
		
		int indexOf_domA = getIndexOfMaximumValue(this.basePredictor.getValueProbability()[a]);
		int indexOf_domB = getIndexOfMaximumValue(this.basePredictor.getValueProbability()[b]);
		String domA = a_uniqueValues.get(indexOf_domA);
		String domB = b_uniqueValues.get(indexOf_domB);
		
		String nonDomA;
		if (a_uniqueValues.size() == 1)
			nonDomA = null;
		else
			nonDomA = a_uniqueValues.get(indexOf_domA == 0 ? 1:0);
		
		String nonDomB;
		if (b_uniqueValues.size() == 1)
			nonDomB = null;
		else
			nonDomB = b_uniqueValues.get(indexOf_domB == 0 ? 1:0);
		
		double p_domA = this.basePredictor.getValueProbability()[a][indexOf_domA];
		double p_domB = this.basePredictor.getValueProbability()[b][indexOf_domB];
		
		double acc_EM = 0.8;
		
		int[] cardinality = numSources_set(aValues, bValues, domA, nonDomA, domB, nonDomB);
		
		double dp_dB_dA = (1 - p_domA) * p_domB * (1 - p_domB) * 
			(-cardinality[1] + cardinality[3] + cardinality[0] - cardinality[2]);
		double dp_dB_nA = (1 - p_domA) * p_domB * (1 - p_domB) * 
				(cardinality[1] - cardinality[3] - cardinality[0] + cardinality[2]);
		
		double dp_dB = p_domA * dp_dB_dA + (1 - p_domA) * dp_dB_nA;
		
		double dp_nB_dA = (1 - p_domA) * (1 - p_domB) * p_domB * 
				(-cardinality[0] + cardinality[2] + cardinality[1] - cardinality[3]);
		double dp_nB_nA = (1 - p_domA) * p_domB * (1 - p_domB) * 
				(cardinality[0] - cardinality[2] - cardinality[1] + cardinality[3]);
			
		double dp_nB = p_domA * dp_nB_dA + (1 - p_domA) * dp_nB_nA;
		
		double dp_B = acc_EM * dp_dB + (1-acc_EM) * dp_nB;
		return dp_B;
	}
	
	/*
	 * function to compute validation by selecting object based on their deltas
	 */
	public double[] computeForDeltas() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		double[] distancesComputed = new double[this.numObjects + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		int count = 0;
		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		accuPR tempPredictor;
		
		int countValidated = 0;
		int countValidatedIterator = 1;
		
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < 50 * countValidatedIterator &&
					!tuplesNotValidated.isEmpty()) {
				double[] delta_p = new double[tuplesNotValidated.size()];
				for (int i = 0; i < tuplesNotValidated.size(); i++) {
					// say object at i has been validated
					int iobj = tuplesNotValidated.get(i);
					for (int j = 0; j < tuplesNotValidated.size(); j++) {
						int jobj = tuplesNotValidated.get(j);
						delta_p[i] += getDelta_b_dueTo_a(iobj, jobj);
					}
				}
				
				int index = getIndexOfMaximumValue(delta_p);
				indices.add(this.truthTuples.get(truthIndices.indexOf(tuplesNotValidated.get(index))));
				tuplesNotValidated.remove(index);
				countValidated++;	
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			
			if (indices.size() > 100 && indices.size() % ((indices.size()/100) * 100) >= 0)
				System.out.println(indices.size() + " : " + countValidated);
			countValidatedIterator++;
		}
				
		return distancesComputed;
	}
	
	/*
	 * args[0] : dataFile, args[1] : truthFile
	 */
	public static void main (String[] args) throws IOException {
		resolvingConflict r = new resolvingConflict(args[0], args[1]);

		int numberOfRuns_s = 1; // for random experiment, how many times to simulate
		double[] distancesComputed_r = r.computeForRandomSelection(numberOfRuns_s);
		double[] distancesComputed_u = r.computeForMU();
		double[] distancesComputed_e = r.computeForMEU();
		double[] distancesComputed_v = r.computeForMVO();
		double[] distancesComputed_d = r.computeForDeltas();
		
		for (int i = 0; i < 10; i++) {
			System.out.println(distancesComputed_r[i] + "\t" + distancesComputed_v[i] + 
					"\t" + distancesComputed_d[i] + "\t" + distancesComputed_e[i] + 
					"\t" + distancesComputed_u[i]);
		}
	}
}
