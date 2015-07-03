import java.io.*;
import java.util.*;

public class syntheticData {
	/*
	 * function to write results to a file
	 */
	public  void writeToFile(String[][] results, String fileName) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			
		for ( int i = 0; i < results.length; i++) {
			for(int j = 0; j < results[i].length; j++ ) {
				if (results[i][j] != null)
					writer.write(results[i][j]);
				if (j != results[i].length - 1)
					writer.write("\t");
			}
			writer.write("\n");
		}
		
		writer.close();
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
	 * function to generate synthetic networked data 
	 * basic grid + hubs + random upto density d
	 */
	public void gridHubsDensity(int p, int q, double x, double d, double max, double min, double a) throws IOException {
		Random r = new Random(System.currentTimeMillis());
		String[][] grid = new String[2*p + 1][2*q + 1];
		int n = 0, m = 0;
		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[i].length; j++) {
				if ((i + j) % 2 != 0)
					grid[i][j] = Integer.toString(m++);
				else {
					if (i % 2 == 0)
						grid[i][j] = Integer.toString(n++);
				}
			}
		}
		
		System.out.println(n + "\t" + m);
		
		// record if object i and source j are connected
		boolean[][] isConnected = new boolean[n][m];
		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[i].length; j++) {
				if ((i % 2 == 0) && (j % 2 == 0)) {
					if (i - 1 >= 0)
						isConnected[Integer.parseInt(grid[i][j])][Integer.parseInt(grid[i-1][j])] = true;
					if (i + 1 < 2*p + 1)
						isConnected[Integer.parseInt(grid[i][j])][Integer.parseInt(grid[i+1][j])] = true;
					if (j - 1 >= 0)
						isConnected[Integer.parseInt(grid[i][j])][Integer.parseInt(grid[i][j-1])] = true;
					if (j + 1 < 2*q + 1)
						isConnected[Integer.parseInt(grid[i][j])][Integer.parseInt(grid[i][j+1])] = true;
				}
			}
		}
		System.out.println(n + "objects and " + m + "sources");
		
		double[] accuracies = new double[m];
//		Arrays.fill(accuracies, a);
		for (int i = 0; i < m; i++) 
			accuracies[i] = r.nextDouble() * (max - min) + min;
			
		String[][] data = new String[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (isConnected[i][j] && data[i][j] == null) {
					if (r.nextDouble() <= accuracies[j])
						data[i][j] = Integer.toString(0);
					else
						data[i][j] = Integer.toString(1);
				}	
			}
			
			List<String> data_row = new ArrayList<String>(Arrays.asList(data[i]));
			data_row.removeAll(Collections.singleton(null));
			List<String> objectValues = new ArrayList<String>(new HashSet<String>(data_row));
			if (objectValues.size() < 2) {
				String hasValue = objectValues.get(0); String putValue = null;
				if (hasValue.equals("1")) // all values are 1
					putValue = Integer.toString(0);
				else
					putValue = Integer.toString(1);
				int location = r.nextInt(m); 
				boolean placed = false;
				while (!placed) {
					if (data[i][location] == null) {
						data[i][location] = putValue;
						placed = !placed;
					}
					else
						location = r.nextInt(m);
				}
			}
		}
	
		String[][] truth = new String[n][2];
		for (int i = 0; i < n; i++) {
			truth[i][0] = Integer.toString(i);
			truth[i][1] = Integer.toString(0);
		}
		
		int countTotal = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (data[i][j] != null) {
					countTotal++;
				}
			}
		}
		System.out.println("Basic grid density = " + (double)countTotal/(m*n));
		
		// add 2 identical hub objects
//		List<Integer> hubObjects = new ArrayList<Integer>();
//		for (int i = 0; i < 5; i++) {
//			int object = r.nextInt(n);
//			hubObjects.add(object);
//			for (int source = 0; source < m; source++)
//				data[object][source] = (r.nextDouble() > accuracies[source] ? Integer.toString(0) : Integer.toString(1));
//		}
		
		// add x% random hub objects
		int hubs = (int) (x * n); 
		int hubDegree = m/4;
		List<Integer> hubObjects = new ArrayList<Integer>();
		for (int i = 0; i < hubs; i++) {
			int object = r.nextInt(n);
			hubObjects.add(object);
			int addedDegree = 0;
			while (addedDegree < hubDegree) {
				int source = r.nextInt(m);
				if (data[object][source] == null) { 
					data[object][source] = (r.nextDouble() > accuracies[source] ? Integer.toString(0) : Integer.toString(1));
					addedDegree++;
				}
			}
		}
		
		// add i votes to hub surroundings
		for (int i = 0; i < hubObjects.size(); i++) {
			int object = hubObjects.get(i);
			for (int j = 0; j < m; j++) {
				if (data[object][j] != null) {
					int added = 0;
					while (added < i) {
						int voteObject = r.nextInt(n);
						if (!hubObjects.contains(voteObject) && data[voteObject][j] == null) {
							data[voteObject][j] = (r.nextDouble() > accuracies[j] ? Integer.toString(0) : Integer.toString(1));
							countTotal++;
							added++;
						}
					}
				}
			}
		}
		System.out.println("Density after adding " + hubObjects.size() + " hubs = " + (double)countTotal/(m*n));
		
		while (countTotal < d * m * n) {
			int object = r.nextInt(n);
			int source = r.nextInt(m);
			if (!hubObjects.contains(object) && data[object][source] == null) {
				data[object][source] = r.nextDouble() < accuracies[source] ? Integer.toString(0) : Integer.toString(1);
				countTotal++;
			}
		}
		System.out.println("Density after adding all edges = " + (double)countTotal/(m*n));
		
		writeToFile(data, "synthetic_data.txt");
		writeToFile(truth, "synthetic_truth.txt");
	}
	
	/*
	 * function to generate synthetic random data
	 */
	public void randomDensity(int n, int m, double d, double min, double max) throws IOException {
		Random r = new Random(System.currentTimeMillis());
		int source, object, countTotal = 0;
		
		String[][] data = new String[n][m];
		String[][] truth = new String[n][2];
		for (int i = 0; i < n; i++) {
			truth[i][0] = Integer.toString(i);
			truth[i][1] = Integer.toString(0);
			data[i][r.nextInt(m)] = Integer.toString(0);
			countTotal++;
		}
		
		double[] accuracies = new double[m];
		for (int i = 0; i < m; i++)
			accuracies[i] = r.nextDouble() * (max - min) + min;
			
		while (countTotal/(double)(m*n) < d) {
			object = r.nextInt(n);
			source = r.nextInt(m);
			if (data[object][source] == null) {
				data[object][source] = (r.nextDouble() < accuracies[source] ? Integer.toString(0) : Integer.toString(1));
				countTotal++;
			}
		}
		
		for (int i = 0; i < n; i++) { // make sure each object has two different values
			List<String> data_row = new ArrayList<String>(Arrays.asList(data[i]));
			data_row.removeAll(Collections.singleton(null));
			List<String> objectValues = new ArrayList<String>(new HashSet<String>(data_row));
			if (objectValues.size() < 2) {
				String hasValue = objectValues.get(0); String putValue = null;
				if (hasValue.equals("1")) // all values are 1
					putValue = Integer.toString(0);
				else if (hasValue.equals("0"))
					putValue = Integer.toString(1);
				int location = r.nextInt(m); 
				boolean placed = false;
				while (!placed) {
					if (data[i][location] == null) {
						data[i][location] = String.valueOf(putValue);
						placed = !placed;
						countTotal++;
					}
					else
						location = r.nextInt(m);
				}
			}
		}
		
		System.out.println(n + " objects \t" + m + " sources \t density: " + (double)countTotal/(m*n));
		writeToFile(data, "synthetic_data.txt");
		writeToFile(truth, "synthetic_truth.txt");
	}
	
	/*
	 * function to generate synthetic chain network
	 */
	public void chain(int n, double d, double max, double min) throws IOException {
		Random r = new Random(System.currentTimeMillis());
		String[][] data = new String[n][n + 1];
		String[][] truth = new String[n][2];
		for (int i = 0; i < n; i++) {
			truth[i][0] = Integer.toString(i);
			truth[i][1] = Integer.toString(0);
		}
		
		double[] accuracies = new double[n + 1];
		for (int i = 0; i < n; i++) 
			accuracies[i] = r.nextDouble() * (max - min) + min;
		
		for (int i = 0; i < n; i++) {
			data[i][i] = String.valueOf(0);
			data[i][i + 1] = String.valueOf(1);
		}
		data[n-1][n] = String.valueOf(1);
		
		int countTotal = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n + 1; j++) {
				if (data[i][j] != null) 
					countTotal++;
			}
		}
		System.out.println(n + " objects \t" + (n + 1) + " sources \t density: " + (double)countTotal/((n + 1)*n));
		
		while (countTotal/(double)((n + 1)*n) < d) {
			int object = r.nextInt(n);
			int source = r.nextInt(n + 1);
			if (data[object][source] == null) {
				data[object][source] = (r.nextDouble() < accuracies[source] ? Integer.toString(0) : Integer.toString(1));
				countTotal++;
			}
		}
		
		System.out.println(n + " objects \t" + (n + 1) + " sources \t density: " + (double)countTotal/((n + 1)*n));
		
		this.writeToFile(data, "chain_data.txt");
		this.writeToFile(truth, "chain_truth.txt");
	}
	
	/*
	 * function to generate network data with hubs for centrality experiments
	 */
	public void centralityGrid(int p, int q, double x) throws IOException {
		// PrintStream out = new PrintStream(new FileOutputStream("output.txt"));
		// System.setOut(out);
		
		Random r = new Random(System.currentTimeMillis());
		int[][] grid = new int[2*p + 1][2*q + 1];
		int n = 0, m = 0; // n objects and m sources
		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[i].length; j++) {
				if ((i + j) % 2 != 0)
					grid[i][j] = m++;
				else {
					if (i % 2 == 0)
						grid[i][j] = n++;
				}
			}
		}
		System.out.print(n + " objects and " + m + " sources\n");
		
		// Generate basic grid: record if object i and source j are connected
		boolean[][] isConnected = new boolean[n][m];
		int countTotal = 0;
		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[i].length; j++) {
				if ((i % 2 == 0) && (j % 2 == 0)) {
					if (i - 1 >= 0)
						isConnected[grid[i][j]][grid[i-1][j]] = true;
					if (i + 1 < 2*p + 1)
						isConnected[grid[i][j]][grid[i+1][j]] = true;
					if (j - 1 >= 0)
						isConnected[grid[i][j]][grid[i][j-1]] = true;
					if (j + 1 < 2*q + 1)
						isConnected[grid[i][j]][grid[i][j+1]] = true;
					countTotal++;
				}
			}
		}
		
		System.out.println("Basic Grid Density = " + countTotal/(double)(m * n));
		
		// inject x% random hub objects
		int hubs = (int) (x * n); 
		int hubDegree = m/3;
		List<Integer> hubObjects = new ArrayList<Integer>();
		for (int i = 0; i < hubs; i++) {
			int object = r.nextInt(n);
			hubObjects.add(object);
			System.out.println("hub " + object);
			int addedDegree = 0;
			while (addedDegree < hubDegree) {
				int source = r.nextInt(m);
				if (!isConnected[object][source]) { 
					isConnected[object][source] = true;
					addedDegree++;
					countTotal++;
				}
			}
		}
		
		double d = countTotal/(double) (m * n);
		System.out.println("Density after adding " + hubs + " hubs = " + d);
		
		while (countTotal/(double)(m*n) < 1.5 * d) {
			int object = r.nextInt(n);
			int source = r.nextInt(m);
			if (!isConnected[object][source] && !hubObjects.contains(object)) {
				isConnected[object][source] = true;
				countTotal++;
			}
		}
		
		for (int i = 0; i < n; i++) {
			int counti = 0;
			for (int j = 0; j < m; j++) 
				if (isConnected[i][j])  
					counti++;
				System.out.println(i + " : " + counti);
		}
		
		System.out.println("Added up to density " + countTotal/(double)(m * n));
		int[] objectCentralities = new int[n]; // #objects reachable by an object
		int[] numSourceVotes = new int[m];		// #votes by a source
		int[] numObjectVotes = new int[n];		// #votes for an object
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m ; j++) { 
				if (isConnected[i][j]) {
					numSourceVotes[j]++;
					numObjectVotes[i]++;
				}
			}
		}

		System.out.println("Network ready");
		
		//----------------------------- network is ready. Now assign votes ------------------------
		
		String[][] data = new String[n][m];
		String one = String.valueOf(1);
		String zero = String.valueOf(0);
		double hubFraction = 0.7; // can tune: probability of having zeroes
		double nonHubFraction = 0.5;
		for (int i = 0; i < n; i++) {
			int countZeroes = 0;
			int zeroes;
			
			if (hubObjects.contains(i))  // assign same entropy to hubs
				zeroes = (int) (hubFraction * numObjectVotes[i]);// can tune;
			else
				zeroes = (int) (nonHubFraction * numObjectVotes[i]);
				
			while (countZeroes < zeroes) {
				int locationOfZero = r.nextInt(m);
				if (isConnected[i][locationOfZero] && data[i][locationOfZero] == null) {
					data[i][locationOfZero] = zero; 
					countZeroes++;
				}
			}
			
			for (int j = 0; j < m; j++)
				if (isConnected[i][j] && data[i][j] == null)
					data[i][j] = one;
							
			// ensure each object has exactly two different values
			List<String> data_row = new ArrayList<String>(Arrays.asList(data[i]));
			data_row.removeAll(Collections.singleton(null));
			List<String> objectValues = new ArrayList<String>(new HashSet<String>(data_row));
			if (objectValues.size() < 2) {
				String hasValue = objectValues.get(0); String putValue = null;
				if (hasValue.equals("1")) // all values are 1
					putValue = zero;
				else
					putValue = one;
				int location = r.nextInt(m); 
				boolean placed = false;
				while (!placed) {
					if (data[i][location] != null) {
						data[i][location] = putValue;
						placed = !placed;
					}
					else
						location = r.nextInt(m);
				}
			}
		}
		
		String[][] truth = new String[n][2];
		for (int i = 0; i < n; i++) {
			truth[i][0] = Integer.toString(i);
			truth[i][1] = Integer.toString(0);
		}
		
		writeToFile(data, "synthetic_data.txt");
		writeToFile(truth, "synthetic_truth.txt");
		
	}

	
	public static void main(String[] args) throws IOException {
		syntheticData g = new syntheticData();
//		g.gridHubsDensity(p, q, x, d, max, min, a);
		int p, q; p = q = 9;
//		g.gridHubsDensity(p, q, 0.02, 0.2, 1, 0.6, 1); // (p x q) grid, x% hubs, d density
		g.randomDensity(300, 10, 0.5, 0.6, 1);
//		g.chain(10, 0.01, 0.6, 1);
//		g.centralityGrid(6, 6, 0.03);
		
		System.out.println("End");
	}
}
