import java.io.*;
import java.util.*;
import java.nio.file.Paths;

public class resolvingConflict {
	public List<List<String>> dataTuples;
	public List<List<String>> truthTuples;
	public int numObjects;
	public int numSources;
	public List<List<String>> objectValues;
	public int[] numberSourceVoted; // array containing the number of objects a source voted for
	public int totalEntries = 0;  // total number of entries in the matrix
	
	public resolvingConflict(String dataFile, String truthFile) throws IOException {
		this.dataTuples   = getTuplesFromFile (dataFile);
		this.truthTuples  = getTuplesFromFile (truthFile);
		
		this.numObjects   = this.dataTuples.size();
		this.numSources   = getNumberOfSources(this.dataTuples);
		
		this.numberSourceVoted = new int[this.numSources];
		getObjectUniqueValues();
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
	            		this.totalEntries++;
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
	 * function to sort an array
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
	 */
	public double computeDatabaseUncertainty(double[][] probabilityValues) {
		double dbUncertainty = 0; //dbUncentrainty = Shannon entropy
		for (int i = 0; i < probabilityValues.length; i++) {
			for (int j = 0; j < probabilityValues[i].length; j++) {
				if (probabilityValues[i][j] > 0)
					dbUncertainty += - probabilityValues[i][j] * Math.log(probabilityValues[i][j]);
			}
		}
		return dbUncertainty;
	}
	
	/*
	 * function to compute centralities of all objects
	 */
	public double[] computeCentralities(List<Integer> tuplesNotValidated) { // for each source count #objects it voted for
		double[] allObjectCentralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++) {
			List<String> objectTupleList = this.dataTuples.get(tuplesNotValidated.get(i));
			for (int j = 0; j < objectTupleList.size(); j++) 
				if (this.dataTuples.get(tuplesNotValidated.get(i)).get(j) != null) // source voted
					allObjectCentralities[i] += this.numberSourceVoted[j];
		}
		return allObjectCentralities;
	}
	
	/*
	 * function to remove an object from further consideration once it has been validated
	 */
	public void decrementSourceCount(int index) {
		for (int i = 0; i < this.dataTuples.get(index).size(); i++) {
			if (this.dataTuples.get(index).get(i) != null)
				this.numberSourceVoted[i]--;
		}
	}
	
	/*
	 * function to compute entropies of all objects
	 */
	public double[] computeEntropies() {
		double fraction = 0;
		double countNotNull = 0;
		List<String> currentValues = new ArrayList<String>();
		double[] allObjectEntropies = new double[this.numObjects];
		
		for (int i = 0; i < this.numObjects; i++) {
			List<String> objectTupleList = this.dataTuples.get(i);
			countNotNull = objectTupleList.size() - Collections.frequency(objectTupleList, null);
			
			currentValues = this.objectValues.get(i);
			
			if (currentValues.size() > 1) {
				for (int j = 0; j < currentValues.size(); j++) {
					fraction = (double) Collections.frequency(objectTupleList, currentValues.get(j)) / countNotNull;
					allObjectEntropies[i]  += -1 * fraction * Math.log(fraction);
				}
			}
		}
		return allObjectEntropies;
	}
	
	/*
	 * function to select object that results in maximum database utility gain (using ground truth)
	 * MUG: Maximum Utility Gain
	 */
	public int selectObjectWithMUG(List<Integer> tuplesNotValidated, List<List<String>> indices,
			accuPR	basePredictor) {
		
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
	public int selectObjectWithMEUG(List<Integer> tuplesNotValidated, List<List<String>> indices,
			accuPR	basePredictor) {
		
		double baseDatabaseUncertainty 	= computeDatabaseUncertainty (basePredictor.getValueProbability());
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
				tupleDBUncertainty += basePredictor.getValueProbability()[tuplesNotValidated.get(i)][j] * 
										computeDatabaseUncertainty(tempValidatedPredictor.getValueProbability());
				denominator += basePredictor.getValueProbability()[tuplesNotValidated.get(i)][j];
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
	 * function to compute validation with selecting object that results in maximum database utility gain
	 * using ground truth
	 */
	public double[] computeForMU(accuPR	basePredictor) throws IOException {
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		accuPR tempPredictor;
		
		distancesComputed[count++] = computeDistanceFromTruth(basePredictor);
		
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		while (!tuplesNotValidated.isEmpty()) {
			// TODO: might want to select with tempPredictor
			int index = selectObjectWithMUG(tuplesNotValidated, indices, basePredictor);
			indices.add(this.truthTuples.get(tuplesNotValidated.get(index)));
			
			tempPredictor = new accuPR(this.dataTuples, indices); 
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			tuplesNotValidated.remove(index);
			if (count%100 == 0)
				System.out.println(count);
		}
		return distancesComputed;
	}
	
	/*
	 * function to compute validation with selecting object that results in maximum expected database utility gain
	 */
	public double[]  computeForMEU(accuPR	basePredictor) {
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>>(); // <index, value true>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		accuPR tempPredictor;
		
		distancesComputed[count++] = computeDistanceFromTruth(basePredictor);
		
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		while (!tuplesNotValidated.isEmpty()) {
			int index = selectObjectWithMEUG(tuplesNotValidated, indices, basePredictor);
			
			indices.add(this.truthTuples.get(tuplesNotValidated.get(index)));
			
			tempPredictor = new accuPR(this.dataTuples, indices);
//			System.out.println(index + ":" + tempPredictor.getIterations());
			
			distancesComputed[count++] += computeDistanceFromTruth(tempPredictor);
			tuplesNotValidated.remove(index);
			if (count%100 == 0)
				System.out.println(count);
		}
		return distancesComputed;
	}
	
	/*
	 * function to compute validation with selecting objects randomly
	 */
	public double[] computeForRandomSelection(accuPR basePredictor, int numberOfRuns) {
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		Random random = new Random(System.currentTimeMillis());
		accuPR tempPredictor;
	    
//		for (int i = 0; i < this.numObjects; i++)
//			tuplesNotValidated.add(i);
		
		distancesComputed[count++] = computeDistanceFromTruth(basePredictor);
		
		for (int i = 0; i < this.truthTuples.size(); i++) {
			if (this.truthTuples.get(i).size() > 1)
				tuplesNotValidated.add(i);
		}
		
		for (int i = 0; i < numberOfRuns; i++) {
			count = 1;
			indices.clear();
			
			int countValidated = 0;
			int countValidatedIterator = 0;
			while (!tuplesNotValidated.isEmpty()) {
				while (indices.size() < Math.pow(2, countValidatedIterator) && 
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
					
			for (int j = 0; j < this.numObjects; j++)
				tuplesNotValidated.add(j);
		}
				
		for (int i = 1; i < distancesComputed.length; i++)
			distancesComputed[i] /= numberOfRuns;
		
		return distancesComputed;
	}
	
	/*
	 * function to compute validation with selecting object based on their entropies
	 * aka MajorityVotingOrdering, or MVO
	 */
	public double[] computeForMVO(accuPR basePredictor) {
		double[] allObjectEntropies = computeEntropies();
		int[] entropyRanks = sortAndRank(allObjectEntropies);
		int index = 0;
		accuPR tempPredictor;
		
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		
		distancesComputed[count++] = computeDistanceFromTruth(basePredictor);
		
		/*
		 * this part for complete dataset truth
		 */
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		/*
		 * START: this part for partial dataset truth
		 */
		
//		double fraction = 0;
//		double countNotNull = 0;
//		List<String> currentValues = new ArrayList<String>();
//		double[] entropies = new double[this.truthTuples.size()];
//		
//		for (int i = 0; i < entropies.length; i++) {
//			List<String> objectTupleList = this.dataTuples.get(Integer.parseInt(this.truthTuples.get(i).get(0)));
//			countNotNull = objectTupleList.size() - Collections.frequency(objectTupleList, null);
//			
//			currentValues = this.objectValues.get(Integer.parseInt(this.truthTuples.get(i).get(0)));
//			
//			if (currentValues.size() > 1) {
//				for (int j = 0; j < currentValues.size(); j++) {
//					fraction = (double) Collections.frequency(objectTupleList, currentValues.get(j)) / countNotNull;
//					allObjectEntropies[i]  += -1 * fraction * Math.log(fraction);
//				}
//			}
//		}
//		
//		entropyRanks = sortAndRank(entropies);
//		
//		for (int i = 0; i < this.truthTuples.size(); i++)
//			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
//		
//		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		/*
		 * END: this part for partial dataset truth
		 */
		
		int countValidated = 0;
		int countValidatedIterator = 13;
		while (!tuplesNotValidated.isEmpty()) {
			while (indices.size() < Math.pow(2, countValidatedIterator) && 
					!tuplesNotValidated.isEmpty()) {
				outerloop:
				for (int i = countValidated; i < entropyRanks.length; i++) {
					if (tuplesNotValidated.contains(entropyRanks[i])) {
						index = tuplesNotValidated.indexOf(entropyRanks[i]);
						break outerloop;
					}
				}
				indices.add(this.truthTuples.get(tuplesNotValidated.get(index)));
//			indices.add(this.truthTuples.get(truthIndices.indexOf(tuplesNotValidated.get(index))));
			
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
	 * function to compute validation with selecting object based on their centralities
	 * or, number of objects voted by sources voting for an object
	 * aka, CentralityOrdering
	 */
	public double[] computeForCO(accuPR basePredictor) {
		double[] distancesComputed = new double[this.numObjects + 1]; // store initial distance from truth at index 0
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int count = 0;
		accuPR tempPredictor;
		
		distancesComputed[count++] = computeDistanceFromTruth(basePredictor);
		
		for (int i = 0; i < this.numObjects; i++)
			tuplesNotValidated.add(i);
		
		
		int countValidated = 0;
		int countValidatedIterator = 0;
		while (!tuplesNotValidated.isEmpty()) {
			
			int addObject;
			
			while (indices.size() < Math.pow(2, countValidatedIterator) 
					&& !tuplesNotValidated.isEmpty()) {
				double[] allObjectCentralities = computeCentralities(tuplesNotValidated);
				
				int[] centralityRanks = sortAndRank(allObjectCentralities);
				
				addObject = tuplesNotValidated.get(centralityRanks[0]);
				indices.add(this.truthTuples.get(addObject));
				decrementSourceCount(addObject);
				tuplesNotValidated.remove(centralityRanks[0]);
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
	 * function to get HL_LH distribution from probability values and accuracies after basic EM
	 * without any validation
	 */
	public void get_HL_LH_distribution(accuPR basePredictor, resolvingConflict r) {
		int HL_indicator = 1; // HL -> High accuracy, low probability
		int LH_indicator = -1; // Low accuracy, high probability
		
		double goodSourceAccuracy_threshold = 0.9;
		double badSourceAccuracy_threshold = 0.7;
		int[][] HL_LH_matrix = new int[this.numObjects][this.numSources];
				
		for (int i = 0; i < this.numObjects; i++) {
			List<String> currentObjectValues = this.objectValues.get(i);
			String value = null;
			if (currentObjectValues.size() != 0) {
				if (basePredictor.getValueProbability()[i][0] > basePredictor.getValueProbability()[i][1])
					value = currentObjectValues.get(0);
				else
					value = currentObjectValues.get(1);
				
				for (int j = 0; j < this.dataTuples.get(i).size(); j++) {
					if (this.dataTuples.get(i).get(j) != null) {
						if (!this.dataTuples.get(i).get(j).equals(value) &&
								basePredictor.getSourceAccuracy()[j] >= goodSourceAccuracy_threshold)
							HL_LH_matrix[i][j] = Integer.valueOf(HL_indicator);
						else if (this.dataTuples.get(i).get(j).equals(value) && 
								basePredictor.getSourceAccuracy()[j] < badSourceAccuracy_threshold)
							HL_LH_matrix[i][j] = Integer.valueOf(LH_indicator);		
					}
				}
			}
		}
		
		// get a sense of interesting votes, i.e. high accuracy source voting for low probability value
		// and low accuracy source voting for high probability value

		int countHL = 0;
		int countLH = 0;
		List<Integer> interestingObjects = new ArrayList<Integer>();
		
		boolean interesting = false;
		for (int i = 0; i < r.numObjects; i++) {
			interesting = false;
			for (int j = 0; j < r.numSources; j++) {
				if (HL_LH_matrix[i][j] == 1) {
					countHL++;
					interesting = true;
				}
				else if (HL_LH_matrix[i][j] == -1) {
					countLH++;
					interesting = true;
				}
			}
			if (interesting)
				interestingObjects.add(i);
		}
		
		System.out.println("% countHL: " + (double)countHL/r.totalEntries);
		System.out.println("% countLH: " + (double)countLH/r.totalEntries);
		System.out.println("#interestingObjects: " + interestingObjects.size());
	}
	
	/*
	 * function to get statistics regarding shortest path distances in object-object network
	 */
	public void get_all_pairs_shortest_distance(resolvingConflict r) {
		int[][] objectsGraph = new int[r.numObjects][r.numObjects];
		for (int i = 0; i < r.numObjects; i++) {
			for (int j = 0; j < r.numObjects; j++) {
				sourceLoop:
				for (int k = 0; k < r.numSources; k++) {
					if (r.dataTuples.get(i).size() > k && r.dataTuples.get(j).size() > k) {
						if (r.dataTuples.get(i).get(k) != null && r.dataTuples.get(j).get(k) != null) {
							objectsGraph[i][j] = 1;
							objectsGraph[j][i] = 1;
							break sourceLoop;
						}
					}
				}
			}
		}
		
		// write to .gml file
		try{
			BufferedWriter writer = new BufferedWriter(new FileWriter("/Users/rpradhan/Desktop/graph.txt"));
//			writer.write("graph [\n");
//			for (int i = 0; i < objectsGraph.length; i++) 
//		    	writer.write("\tnode [\n\t\tid " + i + "\n\t\tlabel \"\"" + i + "\"\"\n\t]\n");

			for (int i = 0; i < objectsGraph.length; i++) {
				for (int j = 0; j < i; j++) {
					if (objectsGraph[i][j] == 1) {
						writer.write(i + "\t" + j + "\t" + 1 + "\n");
					}
				}
			}
//			writer.write("]");
			writer.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	
		int[][] allPairsDistance = new int[r.numObjects][r.numObjects];
		for (int i = 0; i < objectsGraph.length; i++) {
			for (int j = 0; j < objectsGraph.length; j++) {
				if (objectsGraph[i][j] == 1) 
					allPairsDistance[i][j] = 1;
				else
					allPairsDistance[i][j] = Integer.MAX_VALUE - 10;
			}
			allPairsDistance[i][i] = 0;
		}
		
		int[] distanceDistribution = new int[3];
		for (int i = 0; i < objectsGraph.length; i++) {
			for (int j = 0; j < objectsGraph.length; j++) {
				if (i != j && allPairsDistance[i][j] > 1) {
					for (int k = 0; k < objectsGraph.length; k++) {
						if (allPairsDistance[i][j] > allPairsDistance[i][k] + allPairsDistance[k][j]) 
							allPairsDistance[i][j] = allPairsDistance[i][k] + allPairsDistance[k][j];
					}
				}
				System.out.println(allPairsDistance[i][j]);
				if (i < j)
					distanceDistribution[allPairsDistance[i][j]]++;
			}
		}
		
		int max = 0;
		for (int i = 0; i < objectsGraph.length; i++) {
			for (int j = 0; j < objectsGraph.length; j++) {
				if (allPairsDistance[i][j] > max)
					max = allPairsDistance[i][j];
			}
		}
	}
	
	/*
	 * args[0] : dataFile, args[1] : truthFile
	 */
	public static void main (String[] args) throws IOException {
		resolvingConflict r = new resolvingConflict(args[0], args[1]);

		long startTime = System.nanoTime();     // 10^(-9) seconds
		
		// basic EM, no validation
		accuPR	basePredictor = new accuPR (r.dataTuples, null);

		// high accuracy sources voting for low probability values, vice versa
//		r.get_HL_LH_distribution(basePredictor, r);
		
		// all-pairs shortest paths for all objects
//		r.get_all_pairs_shortest_distance(r);
		
		// for small scale experiments
		int numberOfRuns_s = 1; // for random experiment, how many times to simulate
//		double[] distancesComputed = r.computeForRandomSelection(basePredictor, numberOfRuns_s);
//		double[] distancesComputed = r.computeForMU(basePredictor);
//		double[] distancesComputed = r.computeForMEU(basePredictor);

		// for large scale experiments
		int numberOfRuns_l = 1; // for random experiment, how many times to simulate
//		double[] distancesComputed = r.computeForRandomSelection(basePredictor, numberOfRuns_l);
		double[] distancesComputed = r.computeForMVO(basePredictor); // changed MVO to validate objects exponentially
//		double[] distancesComputed = r.computeForCO(basePredictor);

		System.out.println("estimatedTime: " + (System.nanoTime() - startTime)*(Math.pow(10, -9)));
		
		String folderName = (Paths.get(args[0])).getParent().toString();
		r.write_1D_ToFile(distancesComputed, folderName + "/results.txt");
	}
}
