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
			
	public static void main(String[] args) throws IOException {
		syntheticData g = new syntheticData();
		
		Random r1 = new Random(System.currentTimeMillis());
		
		// (p x q) grid : generate n objects and m sources
		int p = 9, q = 9;
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
		System.out.println(n + " objects and " + m + " sources");
		
		double[] accuracies = new double[m];
		int c = 0;

		int numberOfRuns = 1;
		
		for (double initialAccuracy = 0.5; initialAccuracy <= 1; initialAccuracy += 0.1) {
			System.out.println("Accuracy = " + initialAccuracy);
			double[][] values = new double[9][2];
			for (int r = 0; r <= numberOfRuns; r++) {
				Arrays.fill(accuracies, initialAccuracy);
				String[][] data = new String[n][m];
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						if (isConnected[i][j] && data[i][j] == null) {
							c++;
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
						String hasValue = objectValues.get(0); String putValue;
						if (hasValue.equals(1)) // all values are 1
							putValue = Integer.toString(0);
						else
							putValue = Integer.toString(1);
						int location = r1.nextInt(m); 
						boolean placed = false;
						while (!placed) {
							if (data[i][location] == null) {
								data[i][location] = putValue;
								placed = !placed;
								c++;
							}
							else
								location = r1.nextInt(m);
						}
					}
				}
//				System.out.println("Put random votes according to accuracies");
			
				int countTotal = 0;
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						if (data[i][j] != null)
							countTotal++;
					}
				}

				String[][] truth = new String[n][2];
				for (int i = 0; i < n; i++) {
					truth[i][0] = Integer.toString(i);
					truth[i][1] = Integer.toString(0);
				}
									
			// add up to density d
				for (int d = 1; d <= 9; d += 1) {
				// add p random hubs
//				int hubs = 1 + r1.nextInt(10);
//				int hubs = 5;
//				for (int i = 0; i < hubs; i++) {
//					int object = r1.nextInt(n);
//					for (int j = 0; j < m; j++) {
//						if (data[object][j] == null) {
//							data[object][j] = (r1.nextDouble() < accuracies[j] ? Integer.toString(0) : Integer.toString(1));
//							countTotal++;
//						}
//					}
//				}
				
					// randomly add edges upto density d
					int addEdges = ((int)((d * 0.1) * m * n) - countTotal), countAdd = 0;
					addEdges = addEdges >= 0 ? addEdges : 0;
					while (countAdd < addEdges) {
						int object = r1.nextInt(n);
						int source = r1.nextInt(m);
						if (data[object][source] == null) {
							data[object][source] = (r1.nextDouble() < accuracies[source] ? Integer.toString(0) : Integer.toString(1));
							countAdd++;
						}
					}
					countTotal += addEdges;

					g.writeToFile(data, "synthetic_data.txt");
					g.writeToFile(truth, "synthetic_truth.txt");
				
					networkEffect ne = new networkEffect("synthetic_data.txt", "synthetic_truth.txt");
					double[] t = ne.computeForDeltas();
					values[d-1][0] += t[0];
					values[d-1][1] += t[1];
				}
			}
			
			for (int i = 0; i < values.length; i++) { 
				values[i][0] /= values.length;
				values[i][1] /= values.length;
				System.out.println(values[i][0] + "\t" + values[i][1]);
			}
		}
		System.out.println("End \n");
	}
}
