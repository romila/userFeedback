import java.io.*;
import java.util.*;

public class resolvingConflict {
	public List<List<String>> dataTuples;
	public List<List<String>> truthTuples;
	public int numObjects;
	public int numSources;
	public List<List<String>> objectValues; // List of (List of different values for each object)
	public int[] numberSourceVoted; // array containing the number of objects a source voted for
	public accuPR basePredictor;
	public int stepsize = 1;
	public int stepForDistance = 1;
	public int noOfObjectsToValidate = 20;
	public double acc_EM = 0.8;
	public int numberOfHubs = 20;
	
	public resolvingConflict(String dataFile, String truthFile) throws IOException {
		System.out.println("Class constructor");
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
		System.out.println("accuPR converged");
		System.out.println("distance to ground truth " + computeDistanceFromTruth(this.basePredictor));
	}
	
	/*
	 * function to read tuples from a file
	 */
	public List<List<String>> getTuplesFromFile(String fileName) throws IOException {
		System.out.println("Reading from file");
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
	
	/*
	 * function to write a 2D array of values to a file
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
//		System.out.println(distance + "\t" + baseUncertainty);
		return (distance - baseUncertainty)/baseUncertainty;
	}
	
	/*
	 * function to compute overall database uncertainty 
	 * computed as Shannon entropy, U_d = sum over all objects {E}
	 * where E = sum over all values {-p log(p)}
	 */
	public double computeDatabaseUncertainty(double[][] probabilityValues) {
		double dbUncertainty = 0; 
		for (int i = 0; i < probabilityValues.length; i++)
			for (int j = 0; j < probabilityValues[i].length; j++)
				if (probabilityValues[i][j] > 0)
					dbUncertainty += - probabilityValues[i][j] * Math.log(probabilityValues[i][j]);
		return dbUncertainty;
	}
		
	/*
	 * function to compute validation by selecting object that results in minimum distance to ground truth
	 * using ground truth
	 */
	public double[][] computeForMU() throws IOException {
		System.out.println("MU");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, value true>
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		
		accuPR tempPredictor	= this.basePredictor;
		List<List<String>> tempIndices = new ArrayList<List<String>>();
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] distances 	 = new double[tuplesNotValidated.size()];
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				for (int j = 0; j < indices.size(); j++)
					tempIndices.add(indices.get(j));
				tempIndices.add(this.truthTuples.get(tuplesNotValidated.get(i)));
				accuPR	tempValidatedPredictor	= new accuPR(this.dataTuples, tempIndices);
				distances[i] = - computeDistanceFromTruth(tempValidatedPredictor);
				tempIndices.clear();
			}
			
			int[] ranks = sortAndRank(distances);
			indices.add(this.truthTuples.get(tuplesNotValidated.get(ranks[0])));
			tuplesNotValidated.remove(ranks[0]);
			tempPredictor = new accuPR(this.dataTuples, indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
			}
		}
		
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object that results in maximum expected database utility gain
	 */
	public double[][]  computeForMEU() {
		System.out.println("MEU");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
		
		List<List<String>> tempIndices = new ArrayList<List<String>>();
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] expectedUtility = new double[tuplesNotValidated.size()];
			List<String> indexTuple = new ArrayList<String>();
			
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				double tupleDBUncertainty = 0;
				for (int j = 0; j < this.objectValues.get(tuplesNotValidated.get(i)).size(); j++) {
					for (int k = 0; k < indices.size(); k++)
						tempIndices.add(indices.get(k));
					
					indexTuple.clear();
					indexTuple.add(Integer.toString(tuplesNotValidated.get(i))); 
					indexTuple.add(this.objectValues.get(tuplesNotValidated.get(i)).get(j));
					
					tempIndices.add(indexTuple);
					
					accuPR	tempValidatedPredictor	= new accuPR(this.dataTuples, tempIndices);
					tupleDBUncertainty += probabilities[tuplesNotValidated.get(i)][j] * 
											computeDatabaseUncertainty(tempValidatedPredictor.getValueProbability());
					tempIndices.clear();
				}
				expectedUtility[i] = - tupleDBUncertainty;
			}
			
			int[] ranks = sortAndRank(expectedUtility); // object that results in maximum utility
			indices.add(this.truthTuples.get(tuplesNotValidated.get(ranks[0])));
			tuplesNotValidated.remove(ranks[0]);	
			accuPR tempPredictor = new accuPR(this.dataTuples, indices); 
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
			}
		}
				
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting objects randomly
	 */
	public double[][] computeForRandomSelection(int numberOfRuns) {
		System.out.println("Random");
		List<List<String>> indices = new ArrayList<List<String>> (); // list of <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		Random random = new Random(System.currentTimeMillis());
		accuPR tempPredictor;
		
		List<Integer> allTuples = new ArrayList<Integer> ();
		for (int i = 0; i < this.truthTuples.size(); i++) { 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
			allTuples.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		}
		
		for (int i = 0; i < numberOfRuns; i++) {
			indices.clear();
			List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
			
			while (indices.size() < this.noOfObjectsToValidate) {
					int index = random.nextInt(tuplesNotValidated.size());
				    indices.add(this.truthTuples.get(truthIndices.indexOf(tuplesNotValidated.get(index))));
				    tuplesNotValidated.remove(index);
	
				    tempPredictor = new accuPR(this.dataTuples, indices); 
					if (indices.size() % this.stepForDistance == 0) {
				    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
				    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
				    	System.out.println(distancesComputed[count[0]-1][0] + "\t" + distancesComputed[count[1]-1][1]);
					}
			}
					
			// again add all objects to unvalidated list for next run
			tuplesNotValidated.clear();
			for (int j = 0; j < allTuples.size(); j++)
				tuplesNotValidated.add(allTuples.get(i));
		}
				
		for (int i = 1; i < distancesComputed.length; i++) {
			distancesComputed[i][0] /= numberOfRuns;
			distancesComputed[i][1] /= numberOfRuns;
		}
		
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object based on their entropies
	 * aka MajorityVotingOrdering, or MVO
	 */
	public double[][] computeForMVO() {
		System.out.println("MVO");
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
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		List<Integer> tuplesNotValidated = new ArrayList<Integer> ();
		
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		List<Integer> truthIndices = new ArrayList<Integer>(tuplesNotValidated);
		while (indices.size() < this.noOfObjectsToValidate) {
			int start = indices.size();
			int end = start + this.stepsize;
			for (int i = start; i < end; i++) {
				indices.add(this.truthTuples.get(truthIndices.get(entropyRanks[i])));
				tuplesNotValidated.remove(new Integer(entropyRanks[i]));
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
		    	System.out.println(distancesComputed[count[0]-1][0] + "\t" + distancesComputed[count[1]-1][1]);
			}
		}
				
		return distancesComputed;
	}
	
	/*
	 * function to compute validation by selecting object the system has least confidence about
	 */
	public double[][] computeForObjectMEU() {
		System.out.println("MEU_o");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
		
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] obj_meu = new double[tuplesNotValidated.size()];
			for (int i = 0; i < tuplesNotValidated.size(); i++)
				for (int j = 0; j < probabilities[tuplesNotValidated.get(i)].length; j++)
					if (probabilities[tuplesNotValidated.get(i)][j] > 0)
						obj_meu[i] += - probabilities[tuplesNotValidated.get(i)][j] 
								* Math.log(probabilities[tuplesNotValidated.get(i)][j]);
				
			int[] deltaRanks = sortAndRank(obj_meu);

			List<Integer> toDelete = new ArrayList<Integer>();
			for (int i = 0; i < this.stepsize; i++) {
				indices.add(this.truthTuples.get(tuplesNotValidated.get(deltaRanks[i])));
				toDelete.add(tuplesNotValidated.get(deltaRanks[i]));
			}
			
			for (int i = 0; i < toDelete.size(); i++) 
				tuplesNotValidated.remove(new Integer(toDelete.get(i)));
					
			tempPredictor = new accuPR(this.dataTuples, indices);
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
			}
			
			for (int k = 0; k < tuplesNotValidated.size(); k++)
				probabilities[tuplesNotValidated.get(k)] = 
					Arrays.copyOf(tempPredictor.getValueProbability()[tuplesNotValidated.get(k)], 2);
		}
		
		return distancesComputed;
	}
	
	/*
	 * function to validate object based on new centrality measure
	 */
	public double[][] computeForCentrality() {
		System.out.println("Centrality");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++)
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		
		// rank objects according to 1-hop centrality scores.
		double[] newCentralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++)
			for (int j = 0; j < this.dataTuples.get(i).size(); j++)
				if (this.dataTuples.get(i).get(j) != null)  
					newCentralities[i] += this.numberSourceVoted[j];
		
		int[] newRanks = sortAndRank(newCentralities);
		List<Integer> hubs = new ArrayList<Integer>();
		for (int i = 0; i < this.noOfObjectsToValidate; i++)
			hubs.add(newRanks[i]);
		
		for (int i = 0; i < this.noOfObjectsToValidate; i++)
			System.out.println(newCentralities[i]);
		
		accuPR tempPredictor;
		List<List<String>> indices = new ArrayList<List<String>> ();
		for (int i = 0; i < this.noOfObjectsToValidate; i++) {
			indices.add(this.truthTuples.get(hubs.get(i)));
						
			tempPredictor = new accuPR(this.dataTuples, indices); 
			if (indices.size() % this.stepForDistance == 0) {
		    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
			}
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
			nonDomA = domA.equals("1") ? String.valueOf(0):String.valueOf(1);
		
		String nonDomB;
		if (b_uniqueValues.size() == 1)
			nonDomB = null;
		else
			nonDomB = domB.equals("1") ? String.valueOf(0):String.valueOf(1);
		
		double p_domA = pa[indexOf_domA];
		double p_domB = pb[indexOf_domB];
		
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
	
	/*
	 * function to check if two objects have a common source
	 */
	public boolean haveCommonSource(int iobj, int jobj) {
		boolean common = false;
		for (int i = 0; i < this.numSources; i++) {
			if (this.dataTuples.get(iobj).size() > i &&
					this.dataTuples.get(jobj).size() > i &&
					this.dataTuples.get(iobj).get(i) != null &&
					this.dataTuples.get(jobj).get(i) != null) {
				common = true;
				break;
			}
		}
		return common;
	}
	
	/*
	 * approx-network-MEU
	 */
	public double[][] approxNetworkMEU() {
		System.out.println("approx-Network-MEU");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
			
		// rank objects according to 1-hop centrality scores.
		double[] newCentralities = new double[tuplesNotValidated.size()];
		for (int i = 0; i < tuplesNotValidated.size(); i++)
			for (int j = 0; j < this.dataTuples.get(i).size(); j++)
				if (this.dataTuples.get(i).get(j) != null)  
					newCentralities[i] += this.numberSourceVoted[j];
			
		int[] newRanks = sortAndRank(newCentralities);
		List<Integer> hubs = new ArrayList<Integer>();
		for (int i = 0; i < this.numberOfHubs; i++)
			hubs.add(newRanks[i]);
		
		double[][] tempProbabilities = new double[this.numObjects][2];
		System.out.println("Estimated uncertainty reduction \t Q(E)");
		while (indices.size() < this.noOfObjectsToValidate) {
			if(indices.size() < this.numberOfHubs) {
				indices.add(this.truthTuples.get(hubs.get(0)));
				List<Integer> toDelete = new ArrayList<Integer>();
				toDelete.add(hubs.get(0));
				for (int j = 0; j < toDelete.size(); j++) {
					hubs.remove(new Integer(toDelete.get(j)));
					tuplesNotValidated.remove(new Integer(toDelete.get(j)));
				}
			}
			else {
				double[] db_uncertainties = new double[tuplesNotValidated.size()];
				double[] meu_obj = new double[tuplesNotValidated.size()] ;
				for (int i = 0; i < tuplesNotValidated.size(); i++) {
					tempProbabilities = new double[this.numObjects][2];
					for (int k = 0; k < this.numObjects; k++) 
						tempProbabilities[k] = Arrays.copyOf(probabilities[k], 2);
				
					int iobj = tuplesNotValidated.get(i);
					for (int j = 0; j < tuplesNotValidated.size(); j++) {
						int jobj = tuplesNotValidated.get(j);
						if (iobj != jobj && this.haveCommonSource(iobj, jobj))
							tempProbabilities[jobj] = Arrays.copyOf(getDelta(iobj, jobj, 
								tempProbabilities[iobj], tempProbabilities[jobj]), 2);
					}
					Arrays.fill(tempProbabilities[iobj], 0);
					
					for (int j = 0; j < probabilities[iobj].length; j++)
						if (probabilities[iobj][j] > 0)
							meu_obj[i] += - probabilities[iobj][j] * Math.log(probabilities[iobj][j]);
					
					db_uncertainties[i] +=  this.computeDatabaseUncertainty(tempProbabilities)
							- this.computeDatabaseUncertainty(probabilities);
					db_uncertainties[i] /= this.computeDatabaseUncertainty(probabilities);
					
					db_uncertainties[i] *= -1;
					if (db_uncertainties[i] < 0.001)
						db_uncertainties[i] = 0;
				}
					
				int[] deltaRanks = sortAndRank(db_uncertainties);
				int[] meuRanks = sortAndRank(meu_obj);

				List<Integer> toDelete = new ArrayList<Integer>();
				for (int i = 0; i < this.stepsize; i++) {
					if (db_uncertainties[i] > 0) {
						indices.add(this.truthTuples.get(tuplesNotValidated.get(deltaRanks[i])));
						toDelete.add(tuplesNotValidated.get(deltaRanks[i]));
					}
					else {
						indices.add(this.truthTuples.get(tuplesNotValidated.get(meuRanks[i])));
						toDelete.add(tuplesNotValidated.get(meuRanks[i]));
					}
					
					for (int j = 0; j < this.numSources; j++) 
						if (this.dataTuples.get(tuplesNotValidated.get(i)).size() > j &&
								this.dataTuples.get(tuplesNotValidated.get(i)).get(j) != null)
							this.numberSourceVoted[j]--;
				}
				
				for (int i = 0; i < toDelete.size(); i++) 
					tuplesNotValidated.remove(new Integer(toDelete.get(i)));
			}
			
			tempPredictor = new accuPR(this.dataTuples, indices);
			this.basePredictor.setSourceAccuracy(tempPredictor.getSourceAccuracy());
			
			if (indices.size() % this.stepForDistance == 0) {
//				this.qualityOfEstimates(probabilities, tempPredictor.getValueProbability(), indices);
				distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
		    	System.out.println(distancesComputed[count[0]-1][0] + "\t" + distancesComputed[count[1]-1][1]);
			}
			
			for (int k = 0; k < this.numObjects; k++)
				probabilities[k] = Arrays.copyOf(tempPredictor.getValueProbability()[k], 2);
		}
		return distancesComputed;
	}
	
	/*
	 * function to evaluate quality of estimates
	 */
	public void qualityOfEstimates(double[][] probabilities, double[][] actualProbabilities, 
			List<List<String>> indices) {
		double[][] tempProbabilities = new double[this.numObjects][2];
		for (int j = 0; j < this.numObjects; j++)
			tempProbabilities[j] = Arrays.copyOf(probabilities[j], 2);
		
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		for (int i = 0; i < indices.size(); i++)
			tuplesNotValidated.remove(new Integer(indices.get(i).get(0)));
		
		int iobj = Integer.parseInt(indices.get(indices.size() - 1).get(0));
		for (int j = 0; j < tuplesNotValidated.size(); j++) {
			int jobj = tuplesNotValidated.get(j);
			if (iobj != jobj)
				tempProbabilities[jobj] = Arrays.copyOf(getDelta(iobj, jobj,
						tempProbabilities[iobj], tempProbabilities[jobj]), 2);
		}
		
		double abs_distance = 0;
		for (int i = 0; i < this.numObjects; i++)
			for (int j = 0; j < 2; j++)
				if (i != iobj)
					abs_distance += Math.abs(tempProbabilities[i][0] - actualProbabilities[i][0]) +
					Math.abs(tempProbabilities[i][1] - actualProbabilities[i][1]);
		
		abs_distance /= this.numObjects;
		
		System.out.println((this.computeDatabaseUncertainty(tempProbabilities) - 
				this.computeDatabaseUncertainty(this.basePredictor.getValueProbability()))/
				this.computeDatabaseUncertainty(this.basePredictor.getValueProbability()) +
				"\t" + abs_distance);
	}
	
	/*
	 * function to compute validation by selecting object based on their deltas
	 */
	public double[][] approxMEU() {
		System.out.println("approxMEU");
		List<Integer> tuplesNotValidated = new ArrayList<Integer>();
		for (int i = 0; i < this.truthTuples.size(); i++) 
			tuplesNotValidated.add(Integer.parseInt(this.truthTuples.get(i).get(0)));
		
		int[] count = new int[2];
		double[][] distancesComputed = new double[this.numObjects/this.stepsize + 1][2];
		List<List<String>> indices = new ArrayList<List<String>> (); // <index, true value>
		accuPR tempPredictor;
		
		double[][] probabilities = new double[this.numObjects][2];
		for (int i = 0; i < probabilities.length; i++)
			probabilities[i] = Arrays.copyOf(this.basePredictor.getValueProbability()[i], 2);
		
//		System.out.println("Estimated uncertainty reduction \t Q(E)");
		long t1 = System.currentTimeMillis();
		while (indices.size() < this.noOfObjectsToValidate) {
			double[] db_uncertainties = new double[tuplesNotValidated.size()];
			double[] meu_obj = new double[tuplesNotValidated.size()] ;
			for (int i = 0; i < tuplesNotValidated.size(); i++) {
				double[][] tempProbabilities = new double[this.numObjects][2];
				for (int k = 0; k < this.numObjects; k++) 
					tempProbabilities[k] = Arrays.copyOf(probabilities[k], 2);
			
				int iobj = tuplesNotValidated.get(i);
				for (int j = 0; j < tuplesNotValidated.size(); j++) {
					int jobj = tuplesNotValidated.get(j);
					if (iobj != jobj && this.haveCommonSource(iobj, jobj))
						tempProbabilities[jobj] = Arrays.copyOf(getDelta(iobj, jobj, 
							tempProbabilities[iobj], tempProbabilities[jobj]), 2);
				}
				Arrays.fill(tempProbabilities[iobj], 0);
				
				for (int j = 0; j < probabilities[iobj].length; j++)
					if (probabilities[iobj][j] > 0)
						meu_obj[i] += - probabilities[iobj][j] * Math.log(probabilities[iobj][j]);
				
				db_uncertainties[i] +=  this.computeDatabaseUncertainty(tempProbabilities)
						- this.computeDatabaseUncertainty(probabilities);
				db_uncertainties[i] /= this.computeDatabaseUncertainty(probabilities);
				
				db_uncertainties[i] *= -1;
				if (db_uncertainties[i] < 0.001)
					db_uncertainties[i] = 0;
			}
				
			int[] deltaRanks = sortAndRank(db_uncertainties);
			int[] meuRanks = sortAndRank(meu_obj);

			List<Integer> toDelete = new ArrayList<Integer>();
			for (int i = 0; i < this.stepsize; i++) {
				if (db_uncertainties[i] > 0) {
					indices.add(this.truthTuples.get(tuplesNotValidated.get(deltaRanks[i])));
					toDelete.add(tuplesNotValidated.get(deltaRanks[i]));
				}
				else {
					System.out.println("added from meu_o");
					indices.add(this.truthTuples.get(tuplesNotValidated.get(meuRanks[i])));
					toDelete.add(tuplesNotValidated.get(meuRanks[i]));
				}
				
				for (int j = 0; j < this.numSources; j++) 
					if (this.dataTuples.get(tuplesNotValidated.get(i)).size() > j &&
							this.dataTuples.get(tuplesNotValidated.get(i)).get(j) != null)
						this.numberSourceVoted[j]--;
			}
			long t2 = System.currentTimeMillis();
			System.out.println((t2 - t1)/(double)(1000 * this.noOfObjectsToValidate));
			
			for (int i = 0; i < toDelete.size(); i++) 
				tuplesNotValidated.remove(new Integer(toDelete.get(i)));
		
			tempPredictor = new accuPR(this.dataTuples, indices);
			this.basePredictor.setSourceAccuracy(tempPredictor.getSourceAccuracy());
			
			if (indices.size() % this.stepForDistance == 0) {
//				this.qualityOfEstimates(probabilities, tempPredictor.getValueProbability(), indices);
		    	distancesComputed[count[0]++][0] += computeDistanceFromTruth(tempPredictor);
		    	distancesComputed[count[1]++][1] += computeReductionInUncertainty(tempPredictor);
		    	System.out.println(distancesComputed[count[0]-1][0] + "\t" + distancesComputed[count[1]-1][1]);
			}
			
			for (int k = 0; k < this.numObjects; k++)
				probabilities[k] = Arrays.copyOf(tempPredictor.getValueProbability()[k], 2);
		}
		
		return distancesComputed;
	}
	
	/*
	 * args[0] : dataFile, args[1] : truthFile
	 */
	public static void main (String[] args) throws IOException {
//		PrintStream out = new PrintStream(new FileOutputStream("output.txt"));
//		System.setOut(out);
		
		resolvingConflict r = new resolvingConflict(args[0], args[1]);

		int numberOfRuns_s = 1; // for random experiment, how many times to simulate
		long t1 = System.currentTimeMillis();
		
//		double[][] distancesComputed_r = r.computeForRandomSelection(numberOfRuns_s);
//		long t2 = System.currentTimeMillis();
//		double[][] distancesComputed_u = r.computeForMU();
//		long t3 = System.currentTimeMillis();
//		double[][] distancesComputed_e = r.computeForMEU();
//		long t4 = System.currentTimeMillis();
//		double[][] distancesComputed_v = r.computeForMVO();
//		long t5 = System.currentTimeMillis();
//		double[][] distancesComputed_o = r.computeForObjectMEU();
//		long t6 = System.currentTimeMillis();
		double[][] distancesComputed_cc = r.computeForCentrality();
//		double[][] distancesComputed_a = r.approxMEU();
//		long t7 = System.currentTimeMillis();
//		double[][] distancesComputed_n = r.approxNetworkMEU();
		
		// dummies
		double[][] distancesComputed_u = new double[101][2];
		double[][] distancesComputed_e = new double[101][2];
		double[][] distancesComputed_r = new double[101][2];
		double[][] distancesComputed_v = new double[101][2];
		double[][] distancesComputed_a = new double[101][2];
		double[][] distancesComputed_o = new double[101][2];
//		double[][] distancesComputed_cc = new double[101][2];
		double[][] distancesComputed_n = new double[101][2];
		
		System.out.println("Effectiveness");
		System.out.printf("Random\tMU\tMEU\tMVO\tMEU_o\tCentrality\tapprox-MEU\tnetwork-approx-MEU\n");
		for (int i = 0; i < r.noOfObjectsToValidate/r.stepForDistance; i++) 
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
					distancesComputed_r[i][0], distancesComputed_u[i][0], distancesComputed_e[i][0], 
					distancesComputed_v[i][0], distancesComputed_o[i][0], distancesComputed_cc[i][0], 
					distancesComputed_a[i][0], distancesComputed_n[i][0]);
		
		System.out.println("\nReduction in uncertainty");
		System.out.printf("Random\tMU\tMEU\tMVO\tMEU_o\tCentrality\tapprox-MEU\tnetwork-approx-MEU\n");
		for (int i = 0; i < r.noOfObjectsToValidate/r.stepForDistance; i++) 
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
					distancesComputed_r[i][1], distancesComputed_u[i][1], distancesComputed_e[i][1], 
					distancesComputed_v[i][1], distancesComputed_o[i][1], distancesComputed_cc[i][1], 
					distancesComputed_a[i][1], distancesComputed_n[i][1]);
	
//		System.out.println("Efficiency");
//		System.out.println((t4 - t1)/(double)(1000 * r.noOfObjectsToValidate) + "\t" + 
//				(t4 - t4)/(double)(1000 * r.noOfObjectsToValidate) + "\t" +
//				(t4 - t4)/(double)(1000 * r.noOfObjectsToValidate) + "\t"); 
	}
}