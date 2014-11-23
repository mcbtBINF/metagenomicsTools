package mbqc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import utils.ConfigReader;

public class RemoveRare
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(ConfigReader.getMbqcDir()
				+ File.separator + "otuCounts.txt")));
		
		writer.write("otu\tassignment\tcount\n");
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				ConfigReader.getMbqcDir() + File.separator + "merged_otu_filtered.txt")));
		
		String[] headers = reader.readLine().split("\t");
		
		System.out.println(headers.length);
		
		for( int x=0; x < 10; x++)
			System.out.println(headers[x]);
		
		int numDone =0;
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			long count =0;
			for(int x=1; x < splits.length-1; x++)
				count += Long.parseLong(splits[x].replace(".0", ""));
			
			writer.write(splits[0] + "\t" + splits[splits.length-1] + "\t" +  count + "\n");
			
			numDone++;
			
			if( numDone % 100== 0)
			{
				System.out.println(numDone);
				writer.flush();
			}
				
		}
		
		writer.flush(); writer.close();
		reader.close();
	}
}
