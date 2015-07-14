package ruralVsUrbanSvmLight;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;

import utils.ConfigReader;
import utils.ProcessWrapper;

public class LeaveOneOut
{
	private static void writeALine(String s , BufferedWriter writer, boolean mask ,
				int patientToLeaveOut) throws Exception
	{

		String[] splits = s.split("\t");
		
		if( Integer.parseInt(splits[1])== 1 && splits[4].equals("first_A") ) 
			if( mask || Integer.parseInt(splits[2]) != patientToLeaveOut)
			{
				if( mask)
				{
					writer.write( " 0  ");
				}
				else 
				{
					if(splits[3].equals("rural"))
						writer.write("-1  ");
					else if ( splits[3].equals("urban"))
						writer.write("1  ");
					else throw new Exception("No");
				}
				
				for( int x=5; x < splits.length; x++)
				{
					writer.write( x + ":" + splits[x] + " " );
				}
				
				writer.write("\n");
			}
			
			writer.flush();
		}
	
	public static void main(String[] args) throws Exception
	{
		HashSet<Integer> patientIDs = getPatientIDs();
		
		for( Integer i : patientIDs)
		{
			System.out.println("Trying "  + i);
			runATrial(i);
			System.exit(1);
		}
	}
	
	private static String getLeftOut( int patientIDToLeaveOut ) throws Exception
	{
		String returnString = null;
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(ConfigReader.getChinaDir() 
				+ File.separator + "phylum_taxaAsColumnsLogNorm_WithMetadata.txt")));
	
		reader.readLine();
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			if( Integer.parseInt(splits[1])== 1 && splits[4].equals("first_A")
							&& Integer.parseInt(splits[2]) == patientIDToLeaveOut)
			{
				if( returnString != null)
					throw new Exception("No");
				
				returnString = s;
			}
		}
	
		if(returnString ==  null)
			throw new Exception("No");
		
		return returnString;
	}
	
	public static void runATrial(int patientIDToLeaveOut) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(new File(ConfigReader.getChinaDir() 
					+ File.separator + "phylum_taxaAsColumnsLogNorm_WithMetadata.txt")));
		
		reader.readLine();
		
		File outFile = new File(ConfigReader.getSvmDir() + File.separator + 
				"trainWithNo" + patientIDToLeaveOut+ ".txt");
		
		outFile.delete();
		
		if(outFile.exists())
			throw new Exception("Could not delete " + outFile.getAbsolutePath());
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			writeALine(s, writer, false, patientIDToLeaveOut);
		}
		
		writer.flush();  writer.close();
		
		File trainFile = new File(ConfigReader.getSvmDir() + File.separator + 
						"trainedWithout" + patientIDToLeaveOut + ".txt");
		
		
		trainFile.delete();
		
		if( trainFile.exists())
			throw new Exception("No");
		
		// train with svm_learn.exe -z r trainingSet.txt regressModel
		String[] args = new String[3];
		args[0] = ConfigReader.getSvmDir() + File.separator + "svm_learn.exe";
		args[1] = outFile.getAbsolutePath();
		args[2] = trainFile.getAbsolutePath();
		new ProcessWrapper(args);
		
		if( ! trainFile.exists())
			throw new Exception("No");
		
		String leftOut = getLeftOut(patientIDToLeaveOut);
		
		File classifyFile= new File(ConfigReader.getSvmDir() + File.separator +
								"toClassify" + patientIDToLeaveOut + ".txt");
		
		classifyFile.delete();
		
		if( classifyFile.exists())
			throw new Exception("NO");
		
		BufferedWriter classificationWriter = new BufferedWriter(new FileWriter(classifyFile));
		
		writeALine(leftOut, classificationWriter, true, patientIDToLeaveOut);
		
		classificationWriter.flush();  classificationWriter.close();
		
		File svmOut = new File(ConfigReader.getSvmDir() + File.separator +
								"svmResultFor" + patientIDToLeaveOut + ".txt");
		
		svmOut.delete();
		
		if( svmOut.exists())
			throw new Exception("no");
		
		// classify with svm_classify setToClassify.txt regressModel svmOut.txt
		args = new String[4];
		args[0] = ConfigReader.getSvmDir() + File.separator + "svm_classify.exe";
		args[1] = classifyFile.getAbsolutePath();
		args[2] = trainFile.getAbsolutePath();
		args[3] = svmOut.getAbsolutePath();
		new ProcessWrapper(args);
	}
	
	private static HashSet<Integer> getPatientIDs() throws Exception
	{
		HashSet<Integer> patientIDs =new HashSet<Integer>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(ConfigReader.getChinaDir() 
				+ File.separator + "phylum_taxaAsColumnsLogNorm_WithMetadata.txt")));
	
		reader.readLine();
		
		for(String s = reader.readLine() ; s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			patientIDs.add( Integer.parseInt(splits[2]));
		}
		
		return patientIDs;
	}
}