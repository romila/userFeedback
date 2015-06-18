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
	public void gridHubsDensity() throws IOException {
		Random r1 = new Random(System.currentTimeMillis());
		
		// (p x q) grid : generate n objects and m sources
		int p = 21, q = 21;
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
		double initialAccuracy = 0.8;
		Arrays.fill(accuracies, initialAccuracy);
		double x = 0.02; // x% hubs
		
//		for (int i = 0; i < m; i++) {
//			if (i < 600)
//				accuracies[i] = 0.6;
//			else if (i < 1200)
//				accuracies[i] = 0.75;
//			else
//				accuracies[i] = 0.9;
//		}
		
		String[][] data = new String[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (isConnected[i][j] && data[i][j] == null) {
					if (r1.nextDouble() <= accuracies[j])
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
				int location = r1.nextInt(m); 
				boolean placed = false;
				while (!placed) {
					if (data[i][location] == null) {
						data[i][location] = putValue;
						placed = !placed;
					}
					else
						location = r1.nextInt(m);
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
		
		// add x% random hub objects
		int hubs = (int) (x * n); 
		int hubDegree = m/4;
		List<Integer> hubObjects = new ArrayList<Integer>();
		for (int i = 0; i < hubs; i++) {
			int object = r1.nextInt(n);
			hubObjects.add(object);
			int addedDegree = 0;
			while (addedDegree < hubDegree) {
				int source = r1.nextInt(m);
				if (data[object][source] == null) { 
					data[object][source] = (r1.nextDouble() > accuracies[source] ? Integer.toString(0) : Integer.toString(1));
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
						int voteObject = r1.nextInt(n);
						if (!hubObjects.contains(voteObject) && data[voteObject][j] == null) {
							data[voteObject][j] = (r1.nextDouble() > accuracies[j] ? Integer.toString(0) : Integer.toString(1));
							countTotal++;
							added++;
						}
					}
				}
			}
		}
		System.out.println("Density after adding " + hubObjects.size() + " hubs = " + (double)countTotal/(m*n));
		
		// add upto density edges
		double d = 0.1;
		while (countTotal < d * m * n) {
			int object = r1.nextInt(n);
			int source = r1.nextInt(m);
			if (!hubObjects.contains(object) && data[object][source] == null) {
				data[object][source] = r1.nextDouble() < accuracies[source] ? Integer.toString(0) : Integer.toString(1);
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
	public void randomDensity() throws IOException {
		Random r1 = new Random(System.currentTimeMillis());
		int n = 100, m = 10;
		double d = 0.22, a = 0.8;
		int source, object, countTotal = 0;
		
		String[][] data = new String[n][m];
		String[][] truth = new String[n][2];
		for (int i = 0; i < n; i++) {
			truth[i][0] = Integer.toString(i);
			truth[i][1] = Integer.toString(0);
			data[i][r1.nextInt(m)] = Integer.toString(0);
			countTotal++;
		}
		
		while (countTotal/(double)(m*n) < d) {
			object = r1.nextInt(n);
			source = r1.nextInt(m);
			if (data[object][source] == null) {
				data[object][source] = (r1.nextDouble() < a ? Integer.toString(0) : Integer.toString(1));
				countTotal++;
			}
		}
		
		for (int i = 0; i < n; i++) {
			List<String> data_row = new ArrayList<String>(Arrays.asList(data[i]));
			data_row.removeAll(Collections.singleton(null));
			List<String> objectValues = new ArrayList<String>(new HashSet<String>(data_row));
			if (objectValues.size() < 2) {
				String hasValue = objectValues.get(0); String putValue = null;
				if (hasValue.equals("1")) // all values are 1
					putValue = Integer.toString(0);
				else if (hasValue.equals("0"))
					putValue = Integer.toString(1);
				int location = r1.nextInt(m); 
				boolean placed = false;
				while (!placed) {
					if (data[i][location] == null) {
						data[i][location] = String.valueOf(putValue);
						placed = !placed;
					}
					else
						location = r1.nextInt(m);
				}
			}
		}
		
		countTotal = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (data[i][j] != null) 
					countTotal++;
			}
		}
		
		System.out.println(n + " objects \t" + m + " sources \t density: " + (double)countTotal/(m*n));		
		
		writeToFile(data, "synthetic_data.txt");
		writeToFile(truth, "synthetic_truth.txt");
	}
	
	public static void main(String[] args) throws IOException {
		syntheticData g = new syntheticData();
//		g.gridHubsDensity();
		g.randomDensity();
		
		System.out.println("End");
	}
}
