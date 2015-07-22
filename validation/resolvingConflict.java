import java.io.*;
import java.util.*;

public class resolvingConflict {
	public List<List<String>> dataTuples;
	public List<List<String>> truthTuples;
	public int numObjects;
	public int numSources;
	public int numFalseValues;
	public List<List<String>> objectValues; // List of (List of different values for each object)
	public int[] numberSourceVoted; // array containing the number of objects a source voted for
	public accuPR basePredictor;
	public int stepForDistance = 1;
	public int noOfObjectsToValidate = 1;
	public double acc_EM = 0.8;
	public int numberOfHubs = 2;
	public accuPR updated_EM;
	public double baseUncertainty; 
	
	/**
	 * Class constructor
	 * @param dataFile		file of votes to objects by sources
	 * @param truthFile		file of true values of objects
	 * @throws IOException
	 */
	public resolvingConflict(String dataFile, String truthFile) throws IOException {
		this.dataTuples   = getTuplesFromFile (dataFile);
		this.truthTuples  = getTuplesFromFile (truthFile);
		
		this.numObjects   = this.dataTuples.size();
		List<Integer> numSourcesObjects = new ArrayList<Integer>();
		for (int i = 0; i < this.numObjects; i++)
			numSourcesObjects.add(this.dataTuples.get(i).size());
		this.numSources = Collections.max(numSourcesObjects);
		
		this.numberSourceVoted = new int[this.numSources];
		getObjectUniqueValues();
		System.out.println("Running accuPR..");
		
		this.basePredictor = new accuPR (this.dataTuples, null);
		this.baseUncertainty = this.computeDatabaseUncertainty(this.basePredictor.getValueProbability());
		
		System.out.println("accuPR converged");
		System.out.println(this.computeDistanceToTruth(this.basePredictor.valueProbability));
	}
	
	/**
	 * Reads tuples from a file
	 * @param fileName	name of the file to be read
	 * @return			a list of list of object tuples
	 * @throws IOException
	 */
	public List<List<String>> getTuplesFromFile(String fileName) throws IOException {
//		System.out.println("Reading " + fileName);
		try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
			List<List<String>> fileTuples = new ArrayList<List<String>>();
		    String line = reader.readLine();
	        while (line != null) {
	        	String[] attributes = line.split("\t");
	            List<String> singleTuple = new ArrayList<String>();
	            for (int i = 0; i < attributes.length; i++) {
	            	if (!attributes[i].equals(""))
	            		singleTuple.add(attributes[i]);
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
	
	/**
	 * 
	 * @param objectsGraph	array to be written to file
	 * @param fileName		name of file to write array to
	 * @throws IOException
	 */
	public void write_2D_ToFile(double[][] objectsGraph, String fileName) throws IOException {
		System.out.println("Writing to file");
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
			for (int i = 0; i < objectsGraph.length; i++) {
				for (int j = 0; j < objectsGraph[i].length; j++)
					writer.write(objectsGraph[i][j] + "\t");
				writer.write("\n");
			}
			writer.close();
		}
	}
	
	/**
	 * Get unique object values for a tuple
	 */
	public void getObjectUniqueValues() {
		this.objectValues = new ArrayList<List<String>>();
		List<Integer> numUniqueValuesArray = new ArrayList<Integer>();
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
			numUniqueValuesArray.add(objectUniqueValues.size());
		}
		this.numFalseValues	 = Collections.max(numUniqueValuesArray) - 1; // excluding the true value
	}
	
	/**
	 * Get index of maximum value in a set of values
	 * @param	values	double array of values
	 * @return			index of maximum in values
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
	
	/**
	 * Sort an array in decreasing order.
	 * @param	array	array to be sorted
	 * @return			locations of sorted items in array
	 * 					e.g. [4, 5, 3, 6] returns [3, 1, 0, 2]
	 */
	public int[] sortAndRank(double[] array) {
		int j;
		double tempVal;
		int tempRank;
		
		int[] rank = new int[array.length];

		for (int i = 1; i < array.length; i++) {
			j = i;
			rank[i] = i;

			while (j > 0 && array[j-1] < array[j]) {
				tempVal = array[j-1];
				array[j-1] = array[j];
				array[j] = tempVal;

				tempRank = rank[j-1];
				rank[j-1] = rank[j];
				rank[j] = tempRank;

				j--;
			}
		}
		return rank;
	}
	
	/** 
	 * Compute average distance to ground truth assuming ground truth is known for the complete database
	 * @param	probabilities	probabilities of all values of all objects
	 * @return					distanceToTruth d = average over all objects (1 - probability of true value)
	 * *myAngle: d* for one object = sum over all distances from absolute probability distribution
	 * *i.e., if truth = {1, 0, 0} and p = {0.7, 0.2, 0.1}
	 * d* = |1 - 0.7| + |0 - 0.2| + |0 - 0.1| = 0.6, while
	 * d =  |1 - 0.7| = 0.3  
	 * conceptually, d* covers the entire distance but mathematically, while comparing two objects,
	 * d and d* differ only by a factor of 2.
	 */
	public double computeDistanceToTruth(double[][] probabilities) {
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
					distance += Math.abs(1 - probabilities[currentObject][locationOfTrueValue]);
			}
		}
		distance /= (double) this.truthTuples.size();
		return distance;
	}

	/**
	 * Compute percentage reduction in database uncertainty
	 * @param	probabilities	probabilities of all values of all objects
	 * @return					(1 - uncertainty(probabilities)/uncertainty(baseProbabilities))
	 */
	public double computeReductionInUncertainty(double[][] probabilities) {
		double distance = computeDatabaseUncertainty(probabilities);
		return (distance - this.baseUncertainty)/this.baseUncertainty;
	}
	
	/**
	 * Compute database uncertainty as Shannon entropy
	 * @param	probabilities	probabilities of all values of all objects
	 * @return					Shannon entropy of probabilities
	 */
	public double computeDatabaseUncertainty(double[][] probabilities) {
		double dbUncertainty = 0; 
		for (int i = 0; i < probabilities.length; i++)
			for (int j = 0; j < probabilities[i].length; j++)
				if (probabilities[i][j] > 0)
					dbUncertainty += - probabilities[i][j] * Math.log(probabilities[i][j]);
		return dbUncertainty;
	}
		
	/**
	 * Returns probabilities after validating objects in indices and running truth finder
	 * @param	indices		validated objects
	 * @return				probabilities of all values of all objects 
	 */
	public double[][] updateEM( List<List<String>> indices){
		this.updated_EM = new accuPR(this.dataTuples, indices);
		return (this.updated_EM.getValueProbability());
	}
	
	/**
	 * Select object that results in minimum distance to ground truth
	 * @return		resulting distanceToGroundTruth, reductionInUncertainty
	 */
	public double[][] computeGUB() throws IOException {
		System.out.println("GUB");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		List<List<String>> indices = new ArrayList<List<String>> ();
		List<List<String>> tempIndices = new ArrayList<List<String>>();
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] distances 	 = new double[tuplesNotValidated.size()];
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				tempIndices.clear();
				for (int j = 0; j < indices.size(); j++)
					tempIndices.add(indices.get(j));
				List<String> updatedIndex = this.truthTuples.get(tuplesNotValidated.get(i));
				tempIndices.add(updatedIndex);
				newProbabilities = this.updateEM(tempIndices);
				distances[i] = - computeDistanceToTruth(newProbabilities);
			}
			int[] ranks = sortAndRank(distances);
			indices.add(this.truthTuples.get(tuplesNotValidated.get(ranks[0])));
			
			newProbabilities = this.updateEM(indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
			tuplesNotValidated.remove(ranks[0]);
		}
		return distancesComputed;
	}
	
	/**
	 * Select object that results in maximum entropy utility gain
	 * @return		resulting distanceToGroundTruth, reductionInUncertainty
	 */
	public double[][]  computeMEU() {
		System.out.println("MEU");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		int[] count = new int[2];
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		List<List<String>> indices = new ArrayList<List<String>> ();
		List<List<String>> tempIndices = new ArrayList<List<String>>();
		double[][] oldProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		for (int i = 0; i < this.numObjects; i++)
			System.arraycopy(this.basePredictor.getValueProbability()[i], 0, oldProbabilities[i], 0, this.numFalseValues + 1);
		
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] expectedUtility = new double[tuplesNotValidated.size()];
			double baseDBUncertainty = this.computeDatabaseUncertainty(oldProbabilities);
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				int iobj = tuplesNotValidated.get(i);
				double[][] probabilities_0 = new double[this.numObjects][this.numFalseValues + 1];
				double[][] probabilities_1 = new double[this.numObjects][this.numFalseValues + 1];
				for (int j = 0; j < this.objectValues.get(iobj).size(); j++) {
					tempIndices.clear();
					for (int k = 0; k < indices.size(); k++)
						tempIndices.add(indices.get(k));
					
					List<String> tempIndex = new ArrayList<String>();
					tempIndex.add(Integer.toString(iobj)); 
					tempIndex.add(this.objectValues.get(iobj).get(j));
					
					tempIndices.add(tempIndex);
					
					if (j == 0)
						probabilities_0 = this.updateEM(tempIndices);
					if (j == 1)
						probabilities_1 = this.updateEM(tempIndices);
				}
				
				expectedUtility[i] = baseDBUncertainty - 
						(oldProbabilities[iobj][0] * computeDatabaseUncertainty(probabilities_0) +
						 oldProbabilities[iobj][1] * computeDatabaseUncertainty(probabilities_1));

//				// get statistics
//				double[][] stats = new double[this.numObjects][10];
//				
//				// obj-id
//				stats[i][0] = i;
//				
//				// meu
//				stats[i][1] = expectedUtility[i];
//				
//				// meu_o
//				if (oldProbabilities[iobj][0] > 0 && oldProbabilities[iobj][1] > 0)
//					stats[i][2] = - oldProbabilities[iobj][0] * Math.log(oldProbabilities[iobj][0])
//						  	  - oldProbabilities[iobj][1] * Math.log(oldProbabilities[iobj][1]);
//				
//				// meu_rest
//				stats[i][3] = stats[i][1] - stats[i][2];
//				
//				// # 1-hop objects
//				List<Integer> objNeighbor = new ArrayList<Integer>();
//				for (int j = 0; j < this.numObjects; j++) {
//					for (int k = 0; k < this.numSources; k++) {
//						if (i != j && haveCommonSource(i, j)) {
//							objNeighbor.add(j);
//							break;
//						}
//					}
//				}
//				stats[i][4] = objNeighbor.size();
//				
//				// meu_1
//				double meu_1 = 0;
//				double meu_1_before = 0;
//				double meu_1_after = 0;
//				for (int j = 0; j < objNeighbor.size(); j++) {
//					int jobj = objNeighbor.get(j);
//					if (oldProbabilities[jobj][0] > 0 && oldProbabilities[jobj][1] > 0)
//						meu_1_before += - oldProbabilities[jobj][0] * Math.log(oldProbabilities[jobj][0])
//									- oldProbabilities[jobj][1] * Math.log(oldProbabilities[jobj][1]);
//					
//					if (probabilities_0[jobj][0] > 0 && probabilities_0[jobj][1] > 0 &&
//						probabilities_1[jobj][0] > 0 && probabilities_1[jobj][1] > 0)
//					meu_1_after += oldProbabilities[iobj][0] * (- probabilities_0[jobj][0] * Math.log(probabilities_0[jobj][0])
//									- probabilities_0[jobj][1] * Math.log(probabilities_0[jobj][1])) +
//								   oldProbabilities[iobj][0] * (- probabilities_1[jobj][0] * Math.log(probabilities_1[jobj][0])
//									- probabilities_1[jobj][1] * Math.log(probabilities_1[jobj][1]));
//				}
//				stats[i][5] = meu_1_before - meu_1_after;
//				stats[i][6] = stats[i][3] - stats[i][5];
//				
//				// mvo
//				List<String> objectTupleList = this.dataTuples.get(i);
//				int countNotNull = objectTupleList.size() - Collections.frequency(objectTupleList, null);
//				List<String> currentValues = this.objectValues.get(i);
//				
//				if (currentValues.size() > 1) {
//					for (int j = 0; j < currentValues.size(); j++) {
//						double fraction = (double) Collections.frequency(objectTupleList, currentValues.get(j)) 
//								/ (double) countNotNull;
//						if (fraction > 0)
//							stats[i][7]  += - fraction * Math.log(fraction);
//					}
//				}
//				
//				// approx-meu_rest_1hop
//				double[][] tempOldProbabilities = new double[this.numObjects][this.numFalseValues + 1];
//				for (int j = 0; j < this.numObjects; j++)
//					System.arraycopy(oldProbabilities[j], 0, tempOldProbabilities[j], 0, this.numFalseValues);
//					
//				double baseunc = this.computeDatabaseUncertainty(oldProbabilities);
//				for (int j = 0; j < objNeighbor.size(); j++) {
//					int jobj = objNeighbor.get(j);
//					System.arraycopy(getDelta(iobj, jobj, tempOldProbabilities[iobj], tempOldProbabilities[jobj]), 
//						0, tempOldProbabilities[jobj], 0, this.numFalseValues + 1);
//				}
//				stats[i][8] +=  baseunc - this.computeDatabaseUncertainty(tempOldProbabilities);
//				
//				stats[i][9] = baseDBUncertainty - this.computeDatabaseUncertainty(probabilities_1);
//				
//				// approx-meu_rest_2+hop
////				for (int j = 0; j < objNeighbor.size(); j++) {
////					int jobj = objNeighbor.get(j);
////					for (int k = 0; k < this.numObjects; k++) {
////						if (k != i && this.haveCommonSource(jobj, k)) {
////							
////						}
////					}
////				}
//				
//				System.out.println(stats[i][0] + "\t" + stats[i][1] + "\t" + stats[i][2] + "\t" + 
//						stats[i][3] + "\t" + stats[i][4] + "\t" + stats[i][5] + "\t" + stats[i][6] + "\t" + 
//						stats[i][7] + "\t" + stats[i][8] + "\t" + stats[i][9]);
			}
			
			
			int[] ranks = sortAndRank(expectedUtility); // object that results in maximum utility
			indices.add(this.truthTuples.get(tuplesNotValidated.get(ranks[0])));
			newProbabilities = this.updateEM(indices); 
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
			
			for (int k = 0; k < this.numObjects; k++) 
				System.arraycopy(newProbabilities[k], 0, oldProbabilities[k], 0, this.numFalseValues + 1);
			
			tuplesNotValidated.remove(ranks[0]);
		}
				
		return distancesComputed;
	}
	
	/**
	 * Select object at random
	 * @return		resulting distanceToGroundTruth, reductionInUncertainty
	 */
	public double[][] computeRandom(int numberOfRuns) {
		System.out.println("Random");
		List<List<String>> indices = new ArrayList<List<String>> (); // list of <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		Random random = new Random(System.currentTimeMillis());
		
		List<Integer> allTuples = new ArrayList<Integer> ();
		for (int i = 0; i < this.truthTuples.size(); i++) { 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
			allTuples.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		}
		
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		for (int k = 0; k < numberOfRuns; k++) {
			indices.clear();
			
			while (indices.size() < this.noOfObjectsToValidate) {
					int index = random.nextInt(tuplesNotValidated.size());
					List<String> updatedIndex = this.truthTuples.get(tuplesNotValidated.get(index));
					indices.add(updatedIndex);
				    
					newProbabilities = this.updateEM(indices);
					if (indices.size() % this.stepForDistance == 0) {
				    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
				    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
					}
					tuplesNotValidated.remove(index);
			}
					
			// again add all objects to unvalidated list for next run
			tuplesNotValidated.clear();
			for (int j = 0; j < allTuples.size(); j++)
				tuplesNotValidated.add(allTuples.get(j));
		}
				
		for (int i = 1; i < distancesComputed.length; i++) {
			distancesComputed[i][0] /= (double) numberOfRuns;
			distancesComputed[i][1] /= (double) numberOfRuns;
		}
		
		return distancesComputed;
	}
	
	/**
	 * Select object that has highest entropy
	 * @return		resulting distanceToGroundTruth, reductionInUncertainty
	 */
	public double[][] computeLocalMEU() {
		System.out.println("LocaL MEU");
		double[] objectEntropies = new double[this.truthTuples.size()];
		
		// get entropies of objects in the truthfile
		for (int i = 0; i < this.truthTuples.size(); i++) {
			int currentObject = Integer.parseInt(this.truthTuples.get(i).get(0));
			List<String> objectTupleList = this.dataTuples.get(currentObject);
			int countNotNull = objectTupleList.size() - Collections.frequency(objectTupleList, null);
			
			List<String> currentValues = this.objectValues.get(currentObject);
			
			if (currentValues.size() > 1) {
				for (int j = 0; j < currentValues.size(); j++) {
					double fraction = (double) Collections.frequency(objectTupleList, currentValues.get(j)) / (double) countNotNull;
					if (fraction > 0)
						objectEntropies[i]  += - fraction * Math.log(fraction);
				}
			}
		}
		int[] entropyRanks = sortAndRank(objectEntropies); // sort objects in dec order of entropies	
		
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		List<List<String>> indices = new ArrayList<List<String>> ();
		for (int i = 0; i < this.noOfObjectsToValidate; i++) {
			List<String> updatedIndex = this.truthTuples.get(tuplesNotValidated.get(entropyRanks[i]));
			indices.add(updatedIndex);
			newProbabilities = this.updateEM(indices);
			
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
		}
			
		return distancesComputed;
	}
	
	/**
	 * Get number of sources that vote for 1){dA,dB}, 2){dA,nB}, 3){nA,dB}, 4){nA,nB}} 
	 * @param	aValues		votes of sources for object A
	 * @param	bValues		votes of sources for object B
	 * @return				array of number of sources in each of the set above
	 */
	public double[] sources_in_combination(List<String> aValues, List<String> bValues,
								String domA, String nonDomA,
								String domB, String nonDomB) {
		double[] combination = new double[4]; 
		double fraction = 0;
		
		
		// TODO:------------ TO BE UPDATED WITH WITH UPDATED SOURCE ACCURACIES
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
	
	/**
	 * Get change in probabilities of B due to validation of A
	 * @param	object_a	id of object A
	 * @param	object_b	id of object B
	 * @param	probabilities_a	probabilities of values of object A
	 * @param	probabilities_b	probabilities of values of object B
	 * @return					updated probabilities of object B
	 */
	public double[] getDelta(int object_a, int object_b, double[] probabilities_a, double[] probabilities_b) {
		List<String> aValues = this.dataTuples.get(object_a);
		List<String> bValues = this.dataTuples.get(object_b);
		List<String> a_uniqueValues = this.objectValues.get(object_a);
		List<String> b_uniqueValues = this.objectValues.get(object_b);
		
		if (aValues.size() == 0 || bValues.size() == 0)
			return (new double[] {0, 0});
		
		int indexOf_domA = getIndexOfMaximumValue(probabilities_a);
		int indexOf_domB = getIndexOfMaximumValue(probabilities_b);
		
		String domA = null, domB = null;
		if (indexOf_domA < a_uniqueValues.size())
			domA = a_uniqueValues.get(indexOf_domA);
		if (indexOf_domB < b_uniqueValues.size())
			domB = b_uniqueValues.get(indexOf_domB);
		
		String nonDomA;
		if (a_uniqueValues.size() == 1)
			nonDomA = null;
		else
			nonDomA = domA.equals("1") ? String.valueOf(0):String.valueOf(1);
		
		String nonDomB;
		if (b_uniqueValues.size() == 1)
			nonDomB = null;
		else
			nonDomB = domB.equals("1") ? String.valueOf(0):String.valueOf(1);
		
		double p_domA = probabilities_a[indexOf_domA];
		double p_domB = probabilities_b[indexOf_domB];
		
		double[] s_v_tA = sources_in_combination(aValues, bValues, domA, nonDomA, domB, nonDomB);
		
		double dp_dB_dA = (1 - p_domA) * p_domB * (1 - p_domB) * (-s_v_tA[1] + s_v_tA[3] + s_v_tA[0] - s_v_tA[2]);
		double dp_dB_nA =     (p_domA) * p_domB * (1 - p_domB) * ( s_v_tA[1] - s_v_tA[3] - s_v_tA[0] + s_v_tA[2]);
		
		double dp_nB_dA = (1 - p_domA) * (1 - p_domB) * p_domB * (-s_v_tA[0] + s_v_tA[2] + s_v_tA[1] - s_v_tA[3]);
		double dp_nB_nA =     (p_domA) * (1 - p_domB) * p_domB * ( s_v_tA[0] - s_v_tA[2] - s_v_tA[1] + s_v_tA[3]);
				
		double dp_dB = this.acc_EM * dp_dB_dA + (1 - this.acc_EM) * dp_dB_nA;
		double dp_nB = this.acc_EM * dp_nB_dA + (1 - this.acc_EM) * dp_nB_nA;
		
//		if (Math.abs(dp_dB_dA) != 0 || Math.abs(dp_dB_nA) != 0)
//			System.out.println("Halt");
		
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
//		p_dB = Math.abs(p_dB - p_domB); 
//		p_nB = Math.abs(p_nB - (1 - p_domB));
		
//		if (nonDomB == null)
//			return (new double[] {p_dB, p_nB});
		if (domB.equals("1"))
			return (new double[] {p_dB, p_nB}); // probabilities as 2, 1, 0
		else 
			return (new double[] {p_nB, p_dB});
	}
	
	/**
	 * Check if two objects have a common source
	 * @param	object_a	id of object A
	 * @param	object_b	id of object B
	 * @return				true, if A and B have at least one common souce
	 * 						false, otherwise
	 */
	public boolean haveCommonSource(int object_a, int object_b) {
		boolean common = false;
		for (int i = 0; i < this.numSources; i++) {
			if (this.dataTuples.get(object_a).size() > i &&
				this.dataTuples.get(object_b).size() > i &&
				this.dataTuples.get(object_a).get(i) != null &&
				this.dataTuples.get(object_b).get(i) != null) {
				common = true;
				break;
			}
		}
		return common;
	}
	
	/**
	 * Select object that results in maximum uncertainty reduction due to its validation  
	 * by computing approximate probabilities of all objects 
	 */
	public double[][] approxMEU() {
		System.out.println("approxMEU");
		int[] count = new int[2];
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		double[][] oldProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		for (int i = 0; i < this.numObjects; i++)
			System.arraycopy(this.basePredictor.getValueProbability()[i], 0, oldProbabilities[i], 0, this.numFalseValues + 1);
		
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] db_uncertainties = new double[tuplesNotValidated.size()];
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				double[][] tempProbabilities = new double[this.numObjects][2];
				for (int k = 0; k < this.numObjects; k++)
					System.arraycopy(oldProbabilities[k], 0, tempProbabilities[k], 0, this.numFalseValues + 1);
					
				int iobj = tuplesNotValidated.get(i);
				for (int j = 0; j < tuplesNotValidated.size(); j++) {
					int jobj = tuplesNotValidated.get(j);
					if (iobj != jobj && this.haveCommonSource(iobj, jobj))
						System.arraycopy(getDelta(iobj, jobj, tempProbabilities[iobj], tempProbabilities[jobj]), 
								0, tempProbabilities[jobj], 0, this.numFalseValues + 1);
				}
				Arrays.fill(tempProbabilities[iobj], 0);
				
				db_uncertainties[i] +=  this.computeDatabaseUncertainty(tempProbabilities)
						- this.computeDatabaseUncertainty(oldProbabilities);
				if (this.computeDatabaseUncertainty(oldProbabilities) > 0)
					db_uncertainties[i] /= this.computeDatabaseUncertainty(oldProbabilities);
				
				db_uncertainties[i] *= -1;
				if (db_uncertainties[i] < 0.001)
					db_uncertainties[i] = 0;
			}
				
			int[] deltaRanks = sortAndRank(db_uncertainties);
			
			List<String> updatedIndex;
			if (db_uncertainties[0] > 0) {
					updatedIndex = this.truthTuples.get(tuplesNotValidated.get(deltaRanks[0])); 
					tuplesNotValidated.remove(deltaRanks[0]);
			}
			else{
				double[] meu_obj = new double[tuplesNotValidated.size()] ;
				for (int k = 0; k < tuplesNotValidated.size(); k++) {
					int iobj = tuplesNotValidated.get(k);
					for (int j = 0; j < oldProbabilities[iobj].length; j++)
						if (oldProbabilities[iobj][j] > 0)
							meu_obj[k] += - oldProbabilities[iobj][j] * Math.log(oldProbabilities[iobj][j]);
				}
				int[] meuRanks = sortAndRank(meu_obj);

				System.out.println("added from meu_o");
				updatedIndex = this.truthTuples.get(tuplesNotValidated.get(meuRanks[0]));
				tuplesNotValidated.remove(meuRanks[0]);
			}
			
			indices.add(updatedIndex);
			newProbabilities = this.updateEM(indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
			for (int k = 0; k < this.numObjects; k++) 
				System.arraycopy(newProbabilities[k], 0, oldProbabilities[k], 0, this.numFalseValues + 1);
		}
		return distancesComputed;
	}
	
	/**
	 * Estimate the impact of centrality
	 */
	public void validateOne() {
		int validate = 0;
		int oneHopReachability = 0;
		for (int i = 0; i < this.numObjects; i++) {
			List<Integer> objNeighbor = new ArrayList<Integer>();
			for (int j = 0; j < this.numObjects; j++) {
				for (int k = 0; k < this.numSources; k++) {
					if (i != j &&
						this.dataTuples.get(i).size() > k &&
						this.dataTuples.get(j).size() > k &&
						this.dataTuples.get(i).get(k) != null &&
						this.dataTuples.get(j).get(k) != null) { 
						objNeighbor.add(j);
						break;
					}
				}
			}
//			if (i == validate)
//				System.out.println("Centrality of Hub:\t" + objNeighbor.size()); // centrality is number of objects affected
//			else
//				System.out.println(objNeighbor.size());
			if (objNeighbor.size() > oneHopReachability) {
				validate = i;
				oneHopReachability = objNeighbor.size();
			}	
		}
		
		// Compute base
		System.out.println(computeDistanceToTruth(this.basePredictor.getValueProbability()) + "\t" + computeDatabaseUncertainty(this.basePredictor.getValueProbability()));
		
		// Compute with hub validation
		List<List<String>> indices = new ArrayList<List<String>>();
		indices.add(this.truthTuples.get(validate));
		double[][] newProbabilities = new double[this.numObjects][2];
		newProbabilities = this.updateEM(indices);
		System.out.println(computeDistanceToTruth(newProbabilities) + "\t" + computeDatabaseUncertainty(newProbabilities));
		
		
		// alter dataset such that hub has same entropy but different centrality
		int currentObject = Integer.parseInt(this.truthTuples.get(validate).get(0));
		List<String> objectTupleList = this.dataTuples.get(currentObject);
		int countOne = Collections.frequency(objectTupleList, "1");
		int countZero = Collections.frequency(objectTupleList, "0");
		int removeOne = 0;;
		if (countOne > 2) {
			System.out.println("Can remove votes and alter centrality");
			removeOne = 2;
		}
		int removeZero = (int) (removeOne * countZero/countOne);
		
		List<List<String>> tempTuples = new ArrayList<List<String>>();
		Random r = new Random(System.currentTimeMillis());
		for (int i = 0; i < this.dataTuples.size(); i++) {
			if (i != validate)
				tempTuples.add(this.dataTuples.get(i));
			else {
				String[] validateData = this.dataTuples.get(i).toArray(new String[this.dataTuples.get(i).size()]);
				int iterator0 = 0;
				int iterator1 = 0;
				while (iterator0 < removeZero) {
					int location = r.nextInt(this.dataTuples.get(i).size());
					if (validateData[location] != null && validateData[location].equals("0")) {
						validateData[location] = null;
						iterator0++;
					}
				}
				while (iterator1 < removeOne) { 
					int location = r.nextInt(this.dataTuples.get(i).size());
					if (validateData[location] != null && validateData[location].equals("1")) {
						validateData[location] = null;
						iterator1++;
					}
				}
				tempTuples.add(Arrays.asList(validateData));
			}
		}
		this.dataTuples.clear();
		
		for (int i = 0; i < tempTuples.size(); i++)
			this.dataTuples.add(tempTuples.get(i));
		
		List<Integer> objNeighbor = new ArrayList<Integer>();
		for (int j = 0; j < this.numObjects; j++) {
			for (int k = 0; k < this.numSources; k++) {
				if (validate != j &&
					this.dataTuples.get(validate).size() > k &&
					this.dataTuples.get(j).size() > k &&
					this.dataTuples.get(validate).get(k) != null &&
					this.dataTuples.get(j).get(k) != null) { 
					objNeighbor.add(j);
					break;
				}
			}
		}
		
		indices.clear();
		
		// Compute base again
		newProbabilities = this.updateEM(indices);
		System.out.println(computeDistanceToTruth(newProbabilities) + "\t" + computeDatabaseUncertainty(newProbabilities));
		
		// Compute with hub validation
		indices.add(this.truthTuples.get(validate));
		newProbabilities = new double[this.numObjects][2];
		newProbabilities = this.updateEM(indices);
		System.out.println(computeDistanceToTruth(newProbabilities) + "\t" + computeDatabaseUncertainty(newProbabilities));
		System.out.println(oneHopReachability + "\t" + objNeighbor.size());
	}
	
	/**
	 * main() function 
	 * @param	args[0]:dataFile
	 * @param	args[1]:truthFile
	 */
	public static void main (String[] args) throws IOException {
		PrintStream out = new PrintStream(new FileOutputStream("output.txt"));
		System.setOut(out);
		
		resolvingConflict r = new resolvingConflict(args[0], args[1]);

		long t1 = System.currentTimeMillis();
		int numberOfRuns_s = 1; // for random experiment, how many times to simulate
		double[][] distancesComputed_r = r.computeRandom(numberOfRuns_s);
		long t2 = System.currentTimeMillis();
		double[][] distancesComputed_u = r.computeGUB();
		long t3 = System.currentTimeMillis();
		double[][] distancesComputed_e = r.computeMEU();
		long t4 = System.currentTimeMillis();
		double[][] distancesComputed_v = r.computeLocalMEU();
		long t5 = System.currentTimeMillis();
//		double[][] distancesComputed_o = r.computeForObjectMEU();
		long t6 = System.currentTimeMillis();
//		double[][] distancesComputed_cc = r.computeForCentrality();
		long t7 = System.currentTimeMillis();
		double[][] distancesComputed_a = r.approxMEU();
//		double[][] distancesComputed_n = r.approxNetworkMEU();
		long t8 = System.currentTimeMillis();
		
		// dummies
//		double[][] distancesComputed_u = new double[101][2];
//		double[][] distancesComputed_e = new double[101][2];
//		double[][] distancesComputed_r = new double[101][2];
//		double[][] distancesComputed_v = new double[101][2];
		double[][] distancesComputed_o = new double[101][2];
//		double[][] distancesComputed_a = new double[101][2];
		double[][] distancesComputed_cc = new double[101][2];
		double[][] distancesComputed_n = new double[101][2];
		
		System.out.println("Effectiveness");
		System.out.printf("Random\tMU\tMEU\tMVO\tMEU_o\tCentrality\tapprox-MEU\tnetwork-approx-MEU\n");
		for (int i = 0; i < r.noOfObjectsToValidate/r.stepForDistance + 1; i++) 
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
					distancesComputed_r[i][0], distancesComputed_u[i][0], distancesComputed_e[i][0], 
					distancesComputed_v[i][0], distancesComputed_o[i][0], distancesComputed_cc[i][0], 
					distancesComputed_a[i][0], distancesComputed_n[i][0]);
		
		System.out.println("\nReduction in uncertainty");
		System.out.printf("Random\tMU\tMEU\tMVO\tMEU_o\tCentrality\tapprox-MEU\tnetwork-approx-MEU\n");
		for (int i = 0; i < r.noOfObjectsToValidate/r.stepForDistance + 1; i++) 
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
					distancesComputed_r[i][1], distancesComputed_u[i][1], distancesComputed_e[i][1], 
					distancesComputed_v[i][1], distancesComputed_o[i][1], distancesComputed_cc[i][1], 
					distancesComputed_a[i][1], distancesComputed_n[i][1]);
	
		System.out.println("Efficiency");
		System.out.println((t2 - t1)/(double)(1000 * r.noOfObjectsToValidate) + "\t" + 
				(t3 - t2)/(double)(1000 * r.noOfObjectsToValidate) + "\t" +
				(t4 - t3)/(double)(1000 * r.noOfObjectsToValidate) + "\t" +
				(t5 - t4)/(double)(1000 * r.noOfObjectsToValidate) + "\t" +
				(t6 - t5)/(double)(1000 * r.noOfObjectsToValidate) + "\t" + 
				(t7 - t6)/(double)(1000 * r.noOfObjectsToValidate) + "\t" +
				(t8 - t7)/(double)(1000 * r.noOfObjectsToValidate) + "\t");
	}
}
	/**
	 * EXTRA METHODS
	 *//*
	
	*//**
	 * Select object the truth finder has least confidence in
	 * @return		resulting distanceToGroundTruth, reductionInUncertainty
	 *//*
	public double[][] computeForObjectMEU() {
		System.out.println("MEU_o");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		double[][] oldProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		for (int i = 0; i < this.numObjects; i++)
			System.arraycopy(this.basePredictor.getValueProbability()[i], 0, oldProbabilities[i], 0, this.numFalseValues + 1);
		
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] obj_meu = new double[tuplesNotValidated.size()];
			for (int i = 0; i < tuplesNotValidated.size(); i++)
				for (int j = 0; j < oldProbabilities[tuplesNotValidated.get(i)].length; j++)
					if (oldProbabilities[tuplesNotValidated.get(i)][j] > 0)
						obj_meu[i] += - oldProbabilities[tuplesNotValidated.get(i)][j] 
								* Math.log(oldProbabilities[tuplesNotValidated.get(i)][j]);
				
			int[] deltaRanks = sortAndRank(obj_meu);

			List<String> updatedIndex = this.truthTuples.get(tuplesNotValidated.get(deltaRanks[0]));
			indices.add(updatedIndex);
			tuplesNotValidated.remove(deltaRanks[0]);
					
			newProbabilities = this.updateEM(indices);
			if (indices.size() % this.stepForDistance == 0) {
				distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
			
			for (int k = 0; k < this.numObjects; k++) 
				System.arraycopy(newProbabilities[k], 0, oldProbabilities[k], 0, this.numFalseValues + 1);
		}
		
		return distancesComputed;
	}
	
	*//**
	 * Get #objectsToValidate most central objects according to 1-hop score
	 * @param	tuplesNotValidated	list of objects not validated
	 * @return						list of the hubs 
	 *//*
	public List<Integer> getHubs(List<Integer> tuplesNotValidated) {
		double[] centralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < this.numObjects; i++) {
			List<Integer> objNeighbor = new ArrayList<Integer>();
			for (int j = 0; j < this.numObjects; j++) {
				for (int k = 0; k < this.numSources; k++) {
					if (i != j &&
						this.dataTuples.get(i).size() > k &&
						this.dataTuples.get(j).size() > k &&
						this.dataTuples.get(i).get(k) != null &&
						this.dataTuples.get(j).get(k) != null) { 
						objNeighbor.add(j);
						break;
					}
				}
			}
			centralities[i] = objNeighbor.size(); // centrality is number of objects affected
//			System.out.println("Centrality of " + i + "\t" + centralities[i]);
		}
		
		int[] newRanks = sortAndRank(centralities);
		List<Integer> hubs = new ArrayList<Integer>();
		for (int i = 0; i < this.noOfObjectsToValidate; i++)
			hubs.add(tuplesNotValidated.get(newRanks[i]));
		
		return hubs;
	}
	
	*//**
	 * Select object based on centrality measure
	 * @return		resulting distanceToGroundTruth, reductionInUncertainty		
	 *//*
	public double[][] computeForCentrality() {
		System.out.println("Centrality");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		List<Integer> hubs = getHubs(tuplesNotValidated);
		
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		List<List<String>> indices = new ArrayList<List<String>> ();
		for (int i = 0; i < this.noOfObjectsToValidate; i++) {
			List<String> updatedIndex = this.truthTuples.get(hubs.get(i));
			indices.add(updatedIndex);
			
			newProbabilities = this.updateEM(indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
		}
		return distancesComputed;
	}
	
	*//**
	 * approx-network-MEU
	 *//*
	public double[][] approxNetworkMEU() {
		System.out.println("approx-Network-MEU");
		int[] count = new int[2];
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		double[][] oldProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		double[][] newProbabilities = new double[this.numObjects][this.numFalseValues + 1];
		double[][] distancesComputed = new double[this.noOfObjectsToValidate/this.stepForDistance + 1][2];
		distancesComputed[count[0]++][0] += computeDistanceToTruth(this.basePredictor.getValueProbability());
    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(this.basePredictor.getValueProbability());
		
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		for (int i = 0; i < this.numObjects; i++)
			System.arraycopy(this.basePredictor.getValueProbability()[i], 0, oldProbabilities[i], 0, this.numFalseValues + 1);
		
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] sourceVotes = new double[this.numSources];
			for (int i = 0; i < this.numSources; i++) {
				if (this.numberSourceVoted[i] < 2)
					sourceVotes[i] = Integer.MAX_VALUE;
				else
					sourceVotes[i] = this.numberSourceVoted[i];
			}
			
			int[] sourceRanks = sortAndRank(sourceVotes);
			List<Integer> potentialCandidates = new ArrayList<Integer>();
			
			for (int i = this.numSources - 1; i >= (int) (this.numSources * 0.67); i--) {
				int source = sourceRanks[i];
				for (int j = 0; j < this.numObjects; j++) 
					if (tuplesNotValidated.contains(j) && this.dataTuples.get(j).size() > source &&
							this.dataTuples.get(j).get(source) != null)
						potentialCandidates.add(j);
			}
			potentialCandidates = new ArrayList<Integer>(new HashSet<Integer>(potentialCandidates));
			System.out.println(potentialCandidates.size());
			
			double[] db_uncertainties = new double[potentialCandidates.size()];
			for (int i = 0; i < potentialCandidates.size(); i++) {
				int iobj = potentialCandidates.get(i);
				
				double[][] tempProbabilities = new double[this.numObjects][2];
				for (int k = 0; k < this.numObjects; k++)
					System.arraycopy(oldProbabilities[k], 0, tempProbabilities[k], 0, this.numFalseValues + 1);
				
				for (int j = 0; j < tuplesNotValidated.size(); j++) {
					int jobj = tuplesNotValidated.get(j);
					if (iobj != jobj && this.haveCommonSource(iobj, jobj))
						System.arraycopy(getDelta(iobj, jobj, tempProbabilities[iobj], tempProbabilities[jobj]), 
								0, tempProbabilities[jobj], 0, this.numFalseValues + 1);
				}
				Arrays.fill(tempProbabilities[iobj], 0);
				
				db_uncertainties[i] +=  this.computeDatabaseUncertainty(tempProbabilities) 
						- this.computeDatabaseUncertainty(oldProbabilities); 
				db_uncertainties[i] /= this.computeDatabaseUncertainty(oldProbabilities);
				
				db_uncertainties[i] *= -1;
				if (db_uncertainties[i] < 0.001)
					db_uncertainties[i] = 0;
			}
				
			int[] deltaRanks = sortAndRank(db_uncertainties);
			
			List<String> updatedIndex = new ArrayList<String>();
			if (db_uncertainties[0] > 0) {
				updatedIndex = this.truthTuples.get(potentialCandidates.get(deltaRanks[0])); 
				Integer deleteObject = potentialCandidates.get(deltaRanks[0]);
				potentialCandidates.remove(new Integer(deleteObject));
				tuplesNotValidated.remove(new Integer(deleteObject));
			}
			else {
				double[] meu_obj = new double[tuplesNotValidated.size()] ;
				for (int k = 0; k < tuplesNotValidated.size(); k++) {
					int iobj = tuplesNotValidated.get(k);
					for (int j = 0; j < oldProbabilities[iobj].length; j++)
						if (oldProbabilities[iobj][j] > 0)
							meu_obj[k] += - oldProbabilities[iobj][j] * Math.log(oldProbabilities[iobj][j]);
				}
				int[] meuRanks = sortAndRank(meu_obj);
				System.out.println("added from meu_o");
				updatedIndex = this.truthTuples.get(potentialCandidates.get(deltaRanks[0])); 
				Integer deleteObject = tuplesNotValidated.get(meuRanks[0]);
				tuplesNotValidated.remove(new Integer(deleteObject));
			}
			
			indices.add(updatedIndex);
			newProbabilities = this.updateEM(indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceToTruth(newProbabilities);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(newProbabilities);
			}
			for (int k = 0; k < this.numObjects; k++) 
				System.arraycopy(newProbabilities[k], 0, oldProbabilities[k], 0, this.numFalseValues + 1);
		}
		return distancesComputed;
	}
	
	*//**
	 * Evaluate the quality of estimates
	 * @param	probabilities		approximate probabilities 
	 * @param	actualProbabilities	actual probabilities after validation
	 * @param	indices				list of validated objects alongwith true values
	 *//*
	public void qualityOfEstimates(double[][] probabilities, double[][] actualProbabilities, 
			List<List<String>> indices) {
		double[][] tempProbabilities = new double[this.numObjects][2];
		for (int j = 0; j < this.numObjects; j++)
			System.arraycopy(probabilities[j], 0, tempProbabilities[j], 0, this.numFalseValues + 1);
			
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		for (int i = 0; i < indices.size(); i++)
			tuplesNotValidated.remove(new Integer(indices.get(i).get(0)));
		
		int iobj = Integer.parseInt(indices.get(indices.size() - 1).get(0));
		for (int j = 0; j < tuplesNotValidated.size(); j++) {
			int jobj = tuplesNotValidated.get(j);
			if (iobj != jobj)
				System.arraycopy(getDelta(iobj, jobj, tempProbabilities[iobj], tempProbabilities[jobj]), 
						0, tempProbabilities[jobj], 0, this.numFalseValues + 1);
		}
		
		double abs_distance = 0;
		for (int i = 0; i < this.numObjects; i++)
			for (int j = 0; j < 2; j++)
				if (i != iobj)
					abs_distance += Math.abs(tempProbabilities[i][0] - actualProbabilities[i][0]) +
					Math.abs(tempProbabilities[i][1] - actualProbabilities[i][1]);
		
		abs_distance /= (double) this.numObjects;
		
		System.out.println((this.computeDatabaseUncertainty(tempProbabilities) - 
				this.computeDatabaseUncertainty(this.basePredictor.getValueProbability()))/
				this.computeDatabaseUncertainty(this.basePredictor.getValueProbability()) +
				"\t" + abs_distance);
	}
}*/