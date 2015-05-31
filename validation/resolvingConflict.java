import java.io.*;
import java.util.*;

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
	public int stepsize = 1;
	
	public resolvingConflict(String dataFile, String truthFile) throws IOException {
		this.dataTuples   = getTuplesFromFile (dataFile);
		this.truthTuples  = getTuplesFromFile (truthFile);
		
		this.numObjects   = this.dataTuples.size();
		this.numSources   = getNumberOfSources(this.dataTuples);
		
		this.numberSourceVoted = new int[this.numSources];
		getObjectUniqueValues();
		System.out.println("Running accuPR..");
		this.basePredictor = new accuPR (this.dataTuples, null);
		System.out.println("accuPR converged");
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
	 * function to return the reduction in database uncertainty
	 */
	public double computeReductionInUncertainty(accuPR t) {
		double distance = computeDatabaseUncertainty(t.getValueProbability());
		double baseUncertainty = computeDatabaseUncertainty(this.basePredictor.getValueProbability());
		return Math.abs(distance - baseUncertainty);
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
					dbUncertainty += - probabilityValues[i][j] * Math.log(probabilityValues[i][j]);
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
			List<List<String>> tempIndices = new ArrayList<List<String>>(indices);
			tempIndices.add(this.truthTuples.get(tuplesNotValidated.get(i)));
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
		
		accuPR basePredictor	= new accuPR(this.dataTuples, indices);
		double baseDatabaseUncertainty = computeDatabaseUncertainty(basePredictor.getValueProbability());
//		double baseDatabaseUncertainty 	= computeDatabaseUncertainty (this.basePredictor.getValueProbability());
		double[] expectedUtility 	 = new double[tuplesNotValidated.size()];
		
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
				tupleDBUncertainty += basePredictor.getValueProbability()[tuplesNotValidated.get(i)][j] * 
										computeDatabaseUncertainty(tempValidatedPredictor.getValueProbability());
				denominator += basePredictor.getValueProbability()[tuplesNotValidated.get(i)][j];
			}
			tupleDBUncertainty /= denominator;
			
			expectedUtility[i] = - (tupleDBUncertainty - baseDatabaseUncertainty); 
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
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		int count = 0;
		
		accuPR tempPredictor;
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		
		int countValidated = 0;
		int countValidatedIterator = 1;
		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < this.stepsize * countValidatedIterator &&
					!tuplesNotValidated.isEmpty()) {
		
				int index = selectObjectWithMUG(tuplesNotValidated, indices);
			    indices.add(this.truthTuples.get(truthIndices.indexOf(tuplesNotValidated.get(index))));
				tuplesNotValidated.remove(index);
				countValidated++;	
			}
				
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
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
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		int count = 0;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		accuPR tempPredictor;
		
		int countValidated = 0;
		int countValidatedIterator = 1;
		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < this.stepsize * countValidatedIterator &&
					!tuplesNotValidated.isEmpty()) {
				int index = selectObjectWithMEUG(tuplesNotValidated, indices);
			    indices.add(this.truthTuples.get(truthIndices.indexOf(tuplesNotValidated.get(index))));
				tuplesNotValidated.remove(index);
				countValidated++;	
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
		}
				
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting objects randomly
	 */
	public double[] computeForRandomSelection(int numberOfRuns) {
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // list of <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		Random random = new Random(System.currentTimeMillis());
		accuPR tempPredictor;
	    
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
//		distancesComputed[count++] += computeReductionInUncertainty(this.basePredictor);
		System.out.println(distancesComputed[0]);
		
		for (int i = 0; i < this.truthTuples.size(); i++) {
			if (this.truthTuples.get(i).size() > 1)
				tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		}
		
		for (int i = 0; i < numberOfRuns; i++) {
			count = 1;
			indices.clear();
			
			List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
			
			int countValidated = 0;
			int countValidatedIterator = 1;
//			while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < 5) {
				while (indices.size() < this.stepsize * countValidatedIterator  && 
						!tuplesNotValidated.isEmpty()) {
					int index = random.nextInt(tuplesNotValidated.size());
				    indices.add(this.truthTuples.get(truthIndices.indexOf(tuplesNotValidated.get(index))));
				    tuplesNotValidated.remove(index);
					countValidated++;				
				}
				
				tempPredictor = new accuPR(this.dataTuples, indices); 
				distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
//				distancesComputed[count++] += computeReductionInUncertainty(tempPredictor);
				
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
					if (fraction > 0)
						objectEntropies[i]  += - fraction * Math.log(fraction);
				}
			}
		}
		
		// sort objects in decreasing order of their entropies
		int[] entropyRanks = sortAndRank(objectEntropies);		
		accuPR tempPredictor;
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		int countValidated = 0;
		int countValidatedIterator = 1;
//		while (!tuplesNotValidated.isEmpty()) {
		while (indices.size() < 5) {
//			while (indices.size() < this.stepsize * countValidatedIterator &&
//					!tuplesNotValidated.isEmpty()) {
				int start = countValidated;
				int end = start + this.stepsize;
				for (int i = start; i < end; i++) {
					indices.add(this.truthTuples.get(truthIndices.get(entropyRanks[i])));
					tuplesNotValidated.remove(new Integer(entropyRanks[i]));
					countValidated++;
//				}
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
		}
				
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object the system has least confidence about
	 */
	public double[] computeForObjectMEU() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		
		int count = 0;
		int countValidated = 0;
		int countValidatedIterator = 1;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
//		distancesComputed[count++] += computeReductionInUncertainty(this.basePredictor);
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
		
//		while (!tuplesNotValidated.isEmpty()) {
		while (indices.size() < 5) {
//			while (indices.size() < this.stepsize * countValidatedIterator &&
//					!tuplesNotValidated.isEmpty()) {
				
				double[] obj_meu = new double[tuplesNotValidated.size()];
				for (int i = 0; i < tuplesNotValidated.size(); i++) {
					for (int j = 0; j < probabilities[tuplesNotValidated.get(i)].length; j++) {
						if (probabilities[tuplesNotValidated.get(i)][j] > 0)
							obj_meu[i] += - probabilities[tuplesNotValidated.get(i)][j] 
									* Math.log(probabilities[tuplesNotValidated.get(i)][j]);
					}
				}
				
				int[] deltaRanks = sortAndRank(obj_meu);

				List<Integer> toDelete = new ArrayList<Integer>();
				for (int i = 0; i < this.stepsize; i++) {
					indices.add(this.truthTuples.get(tuplesNotValidated.get(deltaRanks[i])));
					toDelete.add(tuplesNotValidated.get(deltaRanks[i]));
					countValidated++;
				}
				
				for (int i = 0; i < toDelete.size(); i++) 
					tuplesNotValidated.remove(new Integer(toDelete.get(i)));
//			}
						
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
//			distancesComputed[count++] += computeReductionInUncertainty(tempPredictor);
			
			for (int k = 0; k < tuplesNotValidated.size(); k++)
				probabilities[tuplesNotValidated.get(k)] = 
					Arrays.copyOf(tempPredictor.getValueProbability()[tuplesNotValidated.get(k)], 2);
		}
		
		return distancesComputed;
	}
	
	/*
	 * function to validate object using centrality method
	 */
	public double[] computeForCentrality() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> ();
		
		int count = 0;
		int countValidated = 0;
		int countValidatedIterator = 1;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
//		distancesComputed[count++] += computeReductionInUncertainty(this.basePredictor);
		
		double[] centralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++) {
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null)
					centralities[i]++;
			}
		}
		
		accuPR tempPredictor;
		int[] deltaRanks = sortAndRank(centralities);
		
		for (int i = 0; i < 10; i++) {
				indices.add(this.truthTuples.get(deltaRanks[i]));
				countValidated++;
						
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
//			distancesComputed[count++] += computeReductionInUncertainty(tempPredictor);
		}
		
		return distancesComputed;
	}
	
	/*
	 * function to validate object based on new centrality measure
	 */
	public double[] computeForNewCentralityScore() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int count = 0;
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		
		double[] centralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++) {
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null) 
					centralities[i]++;	
			}
		}
		
		accuPR tempPredictor;
		int[] deltaRanks = sortAndRank(centralities);
		
		// get hubs - most central nodes
		List<String> hubs = new ArrayList<String>();
		for (int i = 0; i < 10; i++) {
			String object = this.truthTuples.get(deltaRanks[i]).get(0); 
			hubs.add(object);
		}
		
		// rank the hubs according to new scores.
		double[] newCentralities = new double[hubs.size()];
		for (int i = 0; i < hubs.size(); i++) {
			for (int j = 0; j < this.dataTuples.get(Integer.parseInt(hubs.get(i))).size(); j++) {
				if (this.dataTuples.get(Integer.parseInt(hubs.get(i))).get(j) != null)  
					newCentralities[i] += this.numberSourceVoted[j];
			}
		}
		
		int[] newRanks = sortAndRank(newCentralities);
		
		List<List<String>> indices = new ArrayList<List<String>> ();
		for (int i = 0; i < 10; i++) {
			indices.add(this.truthTuples.get(Integer.parseInt(hubs.get(newRanks[i]))));
						
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
		}
		
		return distancesComputed;
	}
	
	/*
	 * for 2-valued case, returns the cardinality of set of sources that vote for
	 * value i in A and j in B, {i,j}
	 * {{dA,dB}, {dA,nB}, {nA,dB}, {nA,nB}}
	 */
	public double[] sources_in_combination(List<String> aValues, List<String> bValues,
								String domA, String nonDomA,
								String domB, String nonDomB) {
		double[] combination = new double[4]; 
		double fraction = 0;
		
		for (int i = 0; i < this.numSources; i++) {
			if (aValues.size() > i && bValues.size() >i &&
					aValues.get(i) != null && bValues.get(i) != null) {
				
				if (this.basePredictor.getSourceAccuracy()[i] != 0 && 
						this.basePredictor.getSourceAccuracy()[i] != 1 &&
						this.numberSourceVoted[i] != 0)
					fraction = 1/(this.numberSourceVoted[i] * this.basePredictor.getSourceAccuracy()[i] * 
						(1 - this.basePredictor.getSourceAccuracy()[i]));
				
				if (aValues.get(i).equals(domA) && bValues.get(i).equals(domB))
					combination[0] += fraction;
				else if (aValues.get(i).equals(domA) && bValues.get(i).equals(nonDomB))
					combination[1] += fraction;
				else if (aValues.get(i).equals(nonDomA) && bValues.get(i).equals(domB))
					combination[2] += fraction;
				else if (aValues.get(i).equals(nonDomA) && bValues.get(i).equals(nonDomB))
					combination[3] += fraction;
			}
		}
		return combination;
	}
	
	/*
	 * function to compute impact on b due to validation of a
	 * 2 - value case
	 */
	public double[] getDelta(int a, int b, double[] pa, double[] pb) {
		List<String> aValues = this.dataTuples.get(a);
		List<String> bValues = this.dataTuples.get(b);
		List<String> a_uniqueValues = this.objectValues.get(a);
		List<String> b_uniqueValues = this.objectValues.get(b);
		
		if (aValues.size() == 0 || bValues.size() == 0)
			return (new double[] {0, 0});
		
		int indexOf_domA = getIndexOfMaximumValue(pa);
		int indexOf_domB = getIndexOfMaximumValue(pb);
		
		String domA = null, domB = null;
		if (indexOf_domA < a_uniqueValues.size())
			domA = a_uniqueValues.get(indexOf_domA);
		if (indexOf_domB < b_uniqueValues.size())
			domB = b_uniqueValues.get(indexOf_domB);
		
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
		
		double p_domA = pa[indexOf_domA];
		double p_domB = pb[indexOf_domB];
		
		double[] s_v_tA = sources_in_combination(aValues, bValues, domA, nonDomA, domB, nonDomB);
		
		double dp_dB_dA = (1 - p_domA) * p_domB * (1 - p_domB) * (-s_v_tA[1] + s_v_tA[3] + s_v_tA[0] - s_v_tA[2]);
		double dp_dB_nA =     (p_domA) * p_domB * (1 - p_domB) * ( s_v_tA[1] - s_v_tA[3] - s_v_tA[0] + s_v_tA[2]);
		
		double dp_nB_dA = (1 - p_domA) * (1 - p_domB) * p_domB * (-s_v_tA[0] + s_v_tA[2] + s_v_tA[1] - s_v_tA[3]);
		double dp_nB_nA =     (p_domA) * (1 - p_domB) * p_domB * ( s_v_tA[0] - s_v_tA[2] - s_v_tA[1] + s_v_tA[3]);
				
		double acc_EM = 0.6;
		double dp_dB = acc_EM * dp_dB_dA + (1 - acc_EM) * dp_dB_nA;
		double dp_nB = acc_EM * dp_nB_dA + (1 - acc_EM) * dp_nB_nA;
		
		// this part is to record the resulting probabilities.
		double p_dB = p_domB + dp_dB;
		if (p_dB > 1)
			p_dB = 1;
		else if (p_dB < 0)
			p_dB = 0;
		
		double p_nB = (1 - p_domB) + dp_nB;
		if (p_nB > 1)
			p_nB = 1;
		else if (p_nB < 0)
			p_nB = 0;
		
		// this part is to record just the change in probabilities. 
//		double p_dB = Math.abs(dp_dB); 
//		double p_nB = Math.abs(dp_nB);
		
		if (nonDomB == null)
			return (new double[] {p_dB, p_nB});
		if (domB.equals("1"))
			return (new double[] {p_dB, p_nB}); // probabilities as 2, 1, 0
		else 
			return (new double[] {p_nB, p_dB});
	}
	
	/*
	 * approx-MEU
	 */
	public double[] approxMEU() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		
		int count = 0;
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
		
		int[][] connected = new int[this.numObjects][this.numObjects];
		for (int i = 0; i < this.numObjects; i++) {
			for (int j = 0; j < this.numObjects; j++) {
				for (int k = 0; k < this.numSources; k++) {
					if (this.dataTuples.get(i).size() > k && this.dataTuples.get(j).size() > k) {
						if (this.dataTuples.get(i).get(k) != null && this.dataTuples.get(j).get(k) != null) {
							connected[i][j] = 1; 
							connected[j][i] = 1;
							break;
						}
					}
				}
			}
		}
		
		while (indices.size() < 5) {
			double[] delta = new double[tuplesNotValidated.size()];
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				int iobj = tuplesNotValidated.get(i);
				for (int j = 0; j < tuplesNotValidated.size(); j++) {
					int jobj = tuplesNotValidated.get(j);
					if (iobj != jobj && connected[iobj][jobj] == 1) {
						double[] netDelta = getDelta(iobj, jobj, probabilities[iobj], probabilities[jobj]);
						delta[i] += netDelta[0] + netDelta[1];
					}
				}
			}
			
			int[] deltaRanks = sortAndRank(delta);
			
			indices.add(this.truthTuples.get(tuplesNotValidated.get(deltaRanks[0])));
			tuplesNotValidated.remove(deltaRanks[0]);
			tempPredictor = new accuPR(this.dataTuples, indices);
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			
			for (int i = 0; i < probabilities.length; i++)
				Arrays.fill(probabilities[i], 0);
			
			for (int k = 0; k < tuplesNotValidated.size(); k++)
				probabilities[tuplesNotValidated.get(k)] = 
					Arrays.copyOf(tempPredictor.getValueProbability()[tuplesNotValidated.get(k)], 2);
		}
		
		return distancesComputed;
	}
	
	/*
	 * approx-network-MEU
	 */
	public double[] approxNetworkMEU() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		
		int count = 0;
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
			
		int[][] connected = new int[this.numObjects][this.numObjects];
		for (int i = 0; i < this.numObjects; i++) {
			for (int j = 0; j < this.numObjects; j++) {
				for (int k = 0; k < this.numSources; k++) {
					if (this.dataTuples.get(i).size() > k && this.dataTuples.get(j).size() > k) {
						if (this.dataTuples.get(i).get(k) != null && this.dataTuples.get(j).get(k) != null) {
							connected[i][j] = 1; 
							connected[j][i] = 1;
							break;
						}
					}
				}
			}
		}
		
		// rank objects according to 1-hop centrality scores.
		double[] newCentralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++) {
			for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
				if (this.dataTuples.get(i).get(j) != null)  
					newCentralities[i] += this.numberSourceVoted[j];
			}
		}
		
		int[] newRanks = sortAndRank(newCentralities);
		List<Integer> hubs = new ArrayList<Integer>();
		for (int i = 0; i < 10; i++)
			hubs.add(newRanks[i]);
		
		while (indices.size() < 5) {
			double[] delta = new double[hubs.size()];
			for (int i = 0; i < delta.length; i++) {
				if (tuplesNotValidated.contains(hubs.get(i))) {
					int iobj = hubs.get(i);
					for (int j = 0; j < tuplesNotValidated.size(); j++) {
						int jobj = tuplesNotValidated.get(j);
						if (iobj != jobj && connected[iobj][jobj] == 1) {
							double[] netDelta = getDelta(iobj, jobj, probabilities[iobj], probabilities[jobj]);
							delta[i] += netDelta[0] + netDelta[1];
						}
					}
				}
//				delta[i] /= newCentralities[i];
			}
			
			int[] deltaRanks = sortAndRank(delta);
			
			indices.add(this.truthTuples.get(hubs.get(deltaRanks[0])));
			tuplesNotValidated.remove(hubs.get(deltaRanks[0]));
			hubs.remove(deltaRanks[0]);
			tempPredictor = new accuPR(this.dataTuples, indices);
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			
			for (int i = 0; i < probabilities.length; i++)
				Arrays.fill(probabilities[i], 0);
			
			for (int k = 0; k < tuplesNotValidated.size(); k++)
				probabilities[tuplesNotValidated.get(k)] = 
					Arrays.copyOf(tempPredictor.getValueProbability()[tuplesNotValidated.get(k)], 2);
		}
		
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object based on their deltas
	 */
	public double[] computeForDeltas() {
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		double[] distancesComputed = new double[this.numObjects/this.stepsize + 1];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		
		int count = 0;
		int countValidated = 0;
		int countValidatedIterator = 1;
		
		distancesComputed[count++] = computeDistanceFromTruth(this.basePredictor);
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
		
		double baseDatabaseUncertainty = computeDatabaseUncertainty(probabilities);
//		System.out.println("--------------------");
		
		while (indices.size() < 5) {
//			while (indices.size() < this.stepsize * countValidatedIterator &&
//					!tuplesNotValidated.isEmpty()) {
//				System.out.println("--------------------");
				
				double[] db_uncertainties = new double[tuplesNotValidated.size()];
				double[] meu_obj = new double[tuplesNotValidated.size()] ;
				for (int i = 0; i < tuplesNotValidated.size(); i++) {
					double[][] tempProbabilities = new double[tuplesNotValidated.size()][2];
					for (int k = 0; k < tuplesNotValidated.size(); k++) 
						tempProbabilities[k] = Arrays.copyOf(probabilities[tuplesNotValidated.get(k)], 2);
				
					int iobj = tuplesNotValidated.get(i);
					for (int j = 0; j < tuplesNotValidated.size(); j++) {
						int jobj = tuplesNotValidated.get(j);
						if (iobj != jobj)
							tempProbabilities[j] = Arrays.copyOf(getDelta(iobj, jobj, 
								tempProbabilities[i], tempProbabilities[j]), 2);
					}
					Arrays.fill(tempProbabilities[i], 0);
					
					for (int j = 0; j < probabilities[iobj].length; j++) {
						if (probabilities[iobj][j] > 0)
							meu_obj[i] += - probabilities[iobj][j] * Math.log(probabilities[iobj][j]);
					}
					
					db_uncertainties[i] =  this.computeDatabaseUncertainty(tempProbabilities)
							- baseDatabaseUncertainty; //final - initial
					
					if (db_uncertainties[i] > 0)
						db_uncertainties[i] = 0;
					
					db_uncertainties[i] *= -1;
//					System.out.println(db_uncertainties[i]);
				}
				
				int[] deltaRanks = sortAndRank(db_uncertainties);

				List<Integer> toDelete = new ArrayList<Integer>();
				for (int i = 0; i < this.stepsize; i++) {
					indices.add(this.truthTuples.get(tuplesNotValidated.get(deltaRanks[i])));
					toDelete.add(tuplesNotValidated.get(deltaRanks[i]));
					countValidated++;
					
					for (int j = 0; j < this.numSources; j++) {
						if (this.dataTuples.get(tuplesNotValidated.get(i)).size() > j &&
								this.dataTuples.get(tuplesNotValidated.get(i)).get(j) != null)
							this.numberSourceVoted[j]--;
					}
				}
				
				for (int i = 0; i < toDelete.size(); i++) 
					tuplesNotValidated.remove(new Integer(toDelete.get(i)));
//			}
		
			tempPredictor = new accuPR(this.dataTuples, indices);
			this.basePredictor.setSourceAccuracy(tempPredictor.getSourceAccuracy());
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			
			for (int i = 0; i < probabilities.length; i++)
				Arrays.fill(probabilities[i], 0);
			
			for (int k = 0; k < tuplesNotValidated.size(); k++)
				probabilities[tuplesNotValidated.get(k)] = 
					Arrays.copyOf(tempPredictor.getValueProbability()[tuplesNotValidated.get(k)], 2);
			baseDatabaseUncertainty = computeDatabaseUncertainty(probabilities);
			
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
		System.out.println("random");
		double[] distancesComputed_r = r.computeForRandomSelection(numberOfRuns_s);
//		double[] distancesComputed_u = r.computeForMU();
//		double[] distancesComputed_e = r.computeForMEU();
		double[] distancesComputed_v = r.computeForMVO();
//		double[] distancesComputed_o = r.computeForObjectMEU();
//		double[] distancesComputed_d = r.computeForDeltas();
//		double[] distancesComputed_c = r.computeForCentrality();
		System.out.println("centrality");
//		double[] distancesComputed_cc = r.computeForNewCentralityScore();
		System.out.println("approx-MEU");
//		double[] distancesComputed_a = r.approxMEU();
		double[] distancesComputed_a = r.computeForDeltas();
		System.out.println("approx-network-MEU");
//		double[] distancesComputed_n = r.approxNetworkMEU();
		
		System.out.println();
		System.out.printf("Random\tCentrality\tapprox-MEU\tapprox-Network-MEU\n");
		for (int i = 0; i < 21; i++) 
			System.out.printf("%.3f\t%.3f\t%.3f\n", distancesComputed_r[i], distancesComputed_v[i], 
					distancesComputed_a[i]);
		
	}
}