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
			
	/* args[0]: n (#objects)
	 * args[1]: m (#sources)
	 * args[2]: a (source accuracy)
	 * args[3]: d (data density) = #entries/mn
	 * Assume: 0 = true value
	 */
	public static void main(String[] args) throws IOException {
		syntheticData g = new syntheticData();
		
		int n = Integer.parseInt(args[0]); n = 300;
		int m = Integer.parseInt(args[1]); m = 5;
		double a = Double.parseDouble(args[2]); a = 0.8;
		double d = Double.parseDouble(args[3]); d = 0.4;
		String[][] data = new String[n][m];
		String[][] truth = new String[n][2];
		
		double[] accuracies = {0.9, 0.85, 0.75, 0.65, 0.6};
		Arrays.fill(accuracies, 0.8);
		
		Random r1 = new Random(System.currentTimeMillis());
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				double one = r1.nextDouble(); // j provided value
				double two = r1.nextDouble(); // j provided true value
				if (one < d) { // j provided value for i
					if (two < accuracies[j]) { // j provided true value for i
						data[i][j] = Integer.toString(0);
					}
					else // source j provided false value for object i
						data[i][j] = Integer.toString(1);
				}
			}			
				
			truth[i][0] = Integer.toString(i);
			truth[i][1] = Integer.toString(0);
		}
		
		g.writeToFile(data, "synthetic_data.txt");
		g.writeToFile(truth, "synthetic_truth.txt");

		
	System.out.println("End \n");
	}
}
