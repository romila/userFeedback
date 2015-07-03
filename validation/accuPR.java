import java.util.*;

public class accuPR {
	public static final double initialAccuracy = 0.8;
	public final double delta 	= 0.001; // error threshold
	public int 		 iterations = 0;
	public List<String> votes 		= new ArrayList<String>();
	public double[][] 	 valueProbability;
	public double[] 	 sourceAccuracy;
	
	public List<String> getVotes(){
		return this.votes;
	}
	
	public double[] getSourceAccuracy(){
		return this.sourceAccuracy;
	}
	
	public void setSourceAccuracy(double[] accuracies){
		for (int i = 0; i < accuracies.length; i++)
			this.sourceAccuracy[i] = accuracies[i];
	}
	
	public double[][] getValueProbability(){
		return this.valueProbability;
	}
	
	public int getIterations(){
		return this.iterations;
	}
	
	/*
	 *  function to round value to N decimal places
	 */
	public static double roundToN(double value, int n){
		return (double)Math.round(value * Math.pow(10, n)) / Math.pow(10, n);
	}
	
	/*
	 *  function to return indices of true values for all objects
	 */
	public static int[] getTrueValueIndices(double[][] a){
		int[] trueValues = new int[a.length];
		double max = 0;
		int index = 0;
	
		for(int i=0; i<a.length; i++){
			max = 0; index = 0;
			for(int j=a[i].length - 1; j>=0; j--){
				if(a[i][j] > max){
					max = a[i][j];
					index = j;
				}
			}
			trueValues[i] = index;				
		}
		return trueValues;
	}

	/* 
	 * function to count number of different elements in two arrays
	 */
	public static int countNumberDifferent(int[] a, int[] b){
		int countDiff = 0;
		for(int i=0; i<a.length; i++){
			if(a[i] != b[i])
				countDiff++;
		}
		return countDiff;
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
	
	public int returnNumberOfSources(List<List<String>> dataTuples){
		List<Integer> numSourcesObjects = new ArrayList<Integer>();
		for(int i=0; i<dataTuples.size(); i++)
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
	
	public accuPR(List<List<String>> dataTuples, List<List<String>> indices) {
		int 				numObjects   = dataTuples.size();
		int 				numSources   = returnNumberOfSources(dataTuples);
		List<List<String>>  objectValues = getObjectUniqueValues(dataTuples);

		List<Integer> numUniqueValuesArray = new ArrayList<Integer>();
        for (int i = 0; i < numObjects; i++)
			numUniqueValuesArray.add(objectValues.get(i).size());
		int 		 numFalseValues		 = Collections.max(numUniqueValuesArray) - 1; // excluding the true value
        
		double[] 	 sourceAccuracyPrev	 = new double[numSources];
		List<String> currentObjectValues = new ArrayList<String>();
		double[]   	 valueConfidence  	 = new double[numFalseValues + 1];
		double 	   	 sumValueProbability = 0;
		int 	   	 countObjects		 = 0;
		List<Integer> truthIndices = new ArrayList<Integer>();
		
		this.valueProbability 	 = new double[numObjects][numFalseValues + 1];
		this.sourceAccuracy		 = new double[numSources];
		
		Arrays.fill(this.sourceAccuracy, initialAccuracy);
		Arrays.fill(sourceAccuracyPrev, -0.1);
			
		if (indices != null) { // some objects have been validated - feed this knowledge into the system
			for (int i = 0; i < indices.size(); i++) {
				int currentObj = Integer.parseInt(indices.get(i).get(0));
				
				truthIndices.add(currentObj);
				currentObjectValues = objectValues.get(currentObj);
				if (indices.get(i).size() > 1) {
					int locationOfTrueValue = currentObjectValues.indexOf(indices.get(i).get(1));
					if (locationOfTrueValue == -1)
						this.valueProbability[currentObj][currentObjectValues.size()] = 1;
					else
						this.valueProbability[currentObj][locationOfTrueValue] = 1;
				}
			}
		}
		
		while (cosineSimilarity(sourceAccuracyPrev, this.sourceAccuracy) < 1 - this.delta) {
			this.iterations++;
			sourceAccuracyPrev = this.sourceAccuracy.clone();
			
			// update value probabilities of objects
			for (int i = 0; i < numObjects; i++) {
				if (!truthIndices.contains(i)) {
					currentObjectValues = objectValues.get(i);
					Arrays.fill(valueConfidence, 0);
					for (int j = 0; j < numSources; j++) 
						if (dataTuples.get(i).size() > j) 
							if (dataTuples.get(i).get(j) != null 
									&& this.sourceAccuracy[j] != 1 
									&& this.sourceAccuracy[j] != 0
									&& numFalseValues > 0) 
								valueConfidence[currentObjectValues.indexOf(dataTuples.get(i).get(j))] += 
										Math.log(numFalseValues * this.sourceAccuracy[j]/(1 - this.sourceAccuracy[j]));
							
					sumValueProbability = 0;
					for (int j = 0; j < numFalseValues + 1; j++) {
						this.valueProbability[i][j] = Math.exp(valueConfidence[j]);
						sumValueProbability += this.valueProbability[i][j];
					}
					
					for (int j = 0; j < numFalseValues + 1; j++)
						this.valueProbability[i][j] =  roundToN(this.valueProbability[i][j]/sumValueProbability, 3);
				}
			}				
			
			// update source accuracies
			Arrays.fill(this.sourceAccuracy, 0);
			for (int i = 0; i < numSources; i++) {
				countObjects = 0;
				for (int j = 0; j < numObjects; j++) {
					currentObjectValues = objectValues.get(j);
					if (dataTuples.get(j).size() > i) {
						if (dataTuples.get(j).get(i) != null) {
							countObjects++;
							this.sourceAccuracy[i] += 
									this.valueProbability[j][currentObjectValues.indexOf(dataTuples.get(j).get(i))];
						}
					}
				}
				this.sourceAccuracy[i] /= (double) countObjects; 
			}
		}
	}
}
