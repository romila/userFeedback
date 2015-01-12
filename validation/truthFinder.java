import java.util.*;

public class truthFinder {
	private final double rho  	= 0.5;
	private final double gamma 	= 0.3;
	private final double delta 	= 0.00001; // error threshold
	private final double t_0 	= 0.9; // initial source accuracy
	private int 		 iterations;
	private List<String> votes;
	private double[][] 	 valueProbability;
	private double[] 	 sourceAccuracy;

	/*
	 * function to get votes for objects
	 */
	public List<String> getVotes() {
		return this.votes;
	}
	
	/*
	 * function to get source accuracies
	 */
	public double[] getSourceAccuracy() {
		return this.sourceAccuracy;
	}
	
	/*
	 * function to get probabilities for different values of objects
	 */
	public double[][] getValueProbability() {
		return this.valueProbability;
	}
	
	/*
	 * function to get number of iterations so far
	 */
	public int getIterations() {
		return this.iterations;
	}
	
	/* 
	 * function to return cosine similarity between two vectors
	 */
	public double cosineSimilarity(double[] a, double[] b) {
		double distance = 0; double sumA = 0; double sumB = 0;
		for (int i = 0; i < a.length; i++) {
			distance += a[i] * b[i]; 
			sumA 	 += a[i] * a[i];
			sumB 	 += b[i] * b[i];
		}
		distance /= Math.sqrt(sumA) * Math.sqrt(sumB);
		return distance;
	}
	
	/*
	 * function to return similarity between two values
	 */
	public double valueSimilarity(double a, double b) {
		return (1/Math.abs(a - b + 1));
	}
	
	/*
	 * function to return number of sources from data tuples
	 */
	public int getNumberOfSources(List<List<String>> dataTuples) {
		List<Integer> numSourcesObjects = new ArrayList<Integer>();
		for (int i = 0; i < dataTuples.size(); i++)
			numSourcesObjects.add(dataTuples.get(i).size());
		return Collections.max(numSourcesObjects);
	}
	
	/* 
	 * function to return unique values for all objects
	 * return objectValues: list of list of values for each object
	 */
	public List<List<String>> getObjectUniqueValues(List<List<String>> dataTuples) {
		List<List<String>> objectValues = new ArrayList<List<String>>();
		for (int i = 0; i < dataTuples.size(); i++) {
			ArrayList<String> objectUniqueValues = new ArrayList<String>();
			for (int j = 0; j < dataTuples.get(i).size(); j++) {
				if(dataTuples.get(i).get(j) != null)
					objectUniqueValues.add(dataTuples.get(i).get(j));
			}
			objectUniqueValues = new ArrayList<String>(new HashSet<String>(objectUniqueValues));
			objectValues.add(objectUniqueValues);
		}
		return objectValues;
	}

	/*
	 * function to return index having maximum in a set of values
	 */
	public int getIndexOfMaximumValue(double[] values) {
		double max = -100; int index = 0; 
		for (int i = 0; i < values.length; i++) {
			if (values[i] > max) {
				max = values[i];
				index = i;
			}
		}
		return index;
	}

	
	/*
	 * function to implement TruthFinder as described in Jiawei Han's work 
	 */
	public truthFinder(List<List<String>> dataTuples, List<List<String>> indices) {
		this.iterations = 0;
		this.votes = new ArrayList<String>();

		int 				numObjects   = dataTuples.size();
		int 				numSources   = getNumberOfSources(dataTuples);
		List<List<String>>  objectValues = getObjectUniqueValues(dataTuples);
		
		List<Integer> numUniqueValuesArray = new ArrayList<Integer>();
        for (int i = 0; i < numObjects; i++)
			numUniqueValuesArray.add(objectValues.get(i).size());
        
       double[] 	 sourceAccuracyPrev	 = new double[numSources];			
		double[] 	 sourceAccuracyTau 	 = new double[numSources];
		double[][] 	 valueConfidence  	 = new double[numObjects][Collections.max(numUniqueValuesArray)];
		List<String> currentObjectValues = new ArrayList<String>();
		
		this.sourceAccuracy	 = new double[numSources];
		this.valueProbability = new double[numObjects][Collections.max(numUniqueValuesArray)];
		
		Arrays.fill(this.sourceAccuracy, this.t_0);
		Arrays.fill(sourceAccuracyPrev, -0.1);
		
		List<Integer> truthIndices = new ArrayList<Integer>();
		
		int[] numberOfObjectsVoted = new int[numSources];
		for (int i = 0; i < numSources; i++) 
			for (int j = 0; j < numObjects; j++) 
				if (dataTuples.get(j).size() > i) 
					if (dataTuples.get(j).get(i) != null)
						numberOfObjectsVoted[i]++;
				
		
		if (indices != null) { // some objects have been validated
			for (int i = 0; i < indices.size(); i++) {
				int currentObj = Integer.parseInt(indices.get(i).get(0));
				truthIndices.add(currentObj);
				currentObjectValues = objectValues.get(currentObj);
				int locationOfTrueValue = currentObjectValues.indexOf(indices.get(i).get(1));
				if (locationOfTrueValue == -1)
					this.valueProbability[currentObj][currentObjectValues.size()] = 1;
				else
					this.valueProbability[currentObj][locationOfTrueValue] = 1;
				
//				for (int j = 0; j < numSources; j++) {
//					if (dataTuples.get(currentObj).size() > j &&
//							dataTuples.get(currentObj).get(j) != null) {
//						if (dataTuples.get(currentObj).get(j).equals(indices.get(i).get(1)))
//							this.sourceAccuracy[j] += (double)1/numberOfObjectsVoted[j];
//						else
//							this.sourceAccuracy[j] -= (double)1/numberOfObjectsVoted[j];
//					}
//				}
			}
		}
		
		// TODO: subtract previous contribution of object, add new contribution
		
		while (cosineSimilarity(sourceAccuracyPrev, this.sourceAccuracy) < 1 - this.delta) {
			this.iterations++;
			sourceAccuracyPrev = this.sourceAccuracy.clone();
			
			for (int i = 0; i < numSources; i++) {
				if (this.sourceAccuracy[i] != 1)
					sourceAccuracyTau[i] = -1 * Math.log(1 - this.sourceAccuracy[i]);
			}

			// update value probabilities of objects
			for (int i = 0; i < numObjects; i++) {
				if (!truthIndices.contains(i)) { // compute only for objects that haven't been validated
					currentObjectValues = objectValues.get(i);
					for (int j = 0; j < numSources; j++) {
						if (dataTuples.get(i).size() > j) {
							if (dataTuples.get(i).get(j) != null)
								valueConfidence[i][currentObjectValues.indexOf(dataTuples.get(i).get(j))] 
										+= sourceAccuracyTau[j];
						}
					}
				}
			}
			
			for (int i = 0; i < numObjects; i++) {
				if (!truthIndices.contains(i)) { 
					for (int j = 0; j < valueConfidence[i].length; j++) {
						double sumWeightedConfidences = 0;
						if(valueConfidence[i][j] != 0) {
							for (int k = 0; k < valueConfidence[i].length; k++) {
								if(j != k)
									sumWeightedConfidences += valueConfidence[i][k] * 
									valueSimilarity(valueConfidence[i][k], valueConfidence[i][j]);
							}
							valueConfidence[i][j] += this.rho * sumWeightedConfidences;
						}
					}
				}
			}
			
			for (int i = 0; i < numObjects; i++) {
				if (!truthIndices.contains(i)) {
					for (int j = 0; j < this.valueProbability[i].length; j++) {
						if (valueConfidence[i][j] != 0)
							this.valueProbability[i][j] = 1/(1 + Math.exp(-1 * this.gamma * valueConfidence[i][j]));
					}
				}
			}

			// update source accuracies
			Arrays.fill(this.sourceAccuracy, 0);
			for (int i = 0; i < numSources; i++) {
				for (int j = 0; j < numObjects; j++) {
					currentObjectValues = objectValues.get(j);
					if (dataTuples.get(j).size() > i) {
						if (dataTuples.get(j).get(i) != null) {
							this.sourceAccuracy[i] += this.valueProbability[j][currentObjectValues.indexOf(dataTuples.get(j).get(i))];
						}
					}
				}
				this.sourceAccuracy[i] /= numberOfObjectsVoted[i]; 
			}
		}

//		 int index = 0;
//		 for (int i = 0; i < numObjects; i++) {
//		 	index = getIndexOfMaximumValue(this.valueProbability[i]);
//			
//		 	if(objectValues.get(i).size() > index)
//		 		this.votes.add(objectValues.get(i).get(index));
//		 	else
//		 		this.votes.add(""); 
//		 }		
	}
}
