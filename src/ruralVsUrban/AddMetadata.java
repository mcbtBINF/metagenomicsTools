package ruralVsUrban;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.StringTokenizer;

import parsers.NewRDPParserFileLine;
import utils.ConfigReader;

public class AddMetadata
{
	private static int getreadNumber(String s)
	{
		return Integer.parseInt(s.split("_")[1].replaceAll("\"", ""));
	}
	
	private static int getPatientId(String s) 
	{
		String first = new StringTokenizer(s, "_").nextToken().replace("A", "").replace("B", "").replace("\"", "");
		return Integer.parseInt(first);
	}
	
	private static String getTimepoint(String s)
	{
		s=s.replaceAll("\"", "");
		if(s.startsWith("B"))
			return "second_B";
		
		return "first_A";
	}
	
	private static void addSomeMetadata(BufferedReader reader, 
					BufferedWriter writer,boolean skipFirst) throws Exception
	{
		writer.write("sampleID\treadNumber\tpatientID\truralUrban\ttimepoint");
		
		if( !skipFirst)
		{
			writer.write("\t" + reader.readLine() + "\n");
		}
		else
		{
			String[] splits = reader.readLine().split("\t");
			for( int x=1; x < splits.length; x++)
				writer.write("\t"  + splits[x]);
			
			writer.write("\n");
		}
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			writer.write(splits[0] + "\t");
			writer.write(getreadNumber(splits[0]) + "\t");
			int patientID = getPatientId(s);
			writer.write(patientID+ "\t");
			
			if( patientID >=1 && patientID <=39)
				writer.write("rural\t");
			else if (patientID >=81 && patientID <= 120)
				writer.write("urban\t");
			else throw new Exception("No");
			
			writer.write(getTimepoint(s));
			
			for( int y=1; y < splits.length; y++)
				writer.write("\t" + splits[y]);
			
			writer.write("\n");	
		}
		
		writer.flush();  writer.close();
		
		reader.close();
	}
	
	/*
	 * Something like:
	 * 
	 *    run the RDP in package cluster stuff on viper
	 *    download the zipped RPP results
	 *    MakeRDPThreeColumnFile
	 *    PivotRDPs
	 *    NormalizeRDPs
	 *    pcoa.txt (in R)
	 *    AddMetadata
	 */
	public static void main(String[] args) throws Exception
	{
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++ )
		{
			System.out.println(NewRDPParserFileLine.TAXA_ARRAY[x]);
			BufferedReader reader = new BufferedReader(new FileReader(new File( 
				ConfigReader.getChinaDir() + File.separator + "pcoa_" + NewRDPParserFileLine.TAXA_ARRAY[x] + ".txt"	)));
			
			BufferedWriter writer =new BufferedWriter(new FileWriter(new File(ConfigReader.getChinaDir() + 
					File.separator + "pcoa_" + NewRDPParserFileLine.TAXA_ARRAY[x] + "_WithMetadata.txt")));
			
			addSomeMetadata(reader, writer,false);
		}
		
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++ )
		{
			System.out.println(NewRDPParserFileLine.TAXA_ARRAY[x]);
			BufferedReader reader = new BufferedReader(new FileReader(new File( 
				ConfigReader.getChinaDir() + File.separator 
				+ NewRDPParserFileLine.TAXA_ARRAY[x] +
				"_taxaAsColumnsLogNorm.txt")));
			
			BufferedWriter writer =new BufferedWriter(new FileWriter(new File(ConfigReader.getChinaDir() + File.separator 
					+ NewRDPParserFileLine.TAXA_ARRAY[x] +
					"_taxaAsColumnsLogNorm" + "_WithMetadata.txt")));
			
			addSomeMetadata(reader, writer,true);
		}
	}
}
