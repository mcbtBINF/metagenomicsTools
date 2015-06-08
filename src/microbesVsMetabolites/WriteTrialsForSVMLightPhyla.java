package microbesVsMetabolites;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;

import microbesVsMetabolites.WriteTrialsForSVMLight.MetaboliteClass;
import parsers.OtuWrapper;
import utils.ConfigReader;
import utils.Pearson;
import utils.ProcessWrapper;
import utils.Regression;
import utils.TabReader;

public class WriteTrialsForSVMLightPhyla
{
	
	private static void deleteOrThrow(File file) throws Exception
	{
		file.delete();
		
		if( file.exists())
			throw new Exception("Could not delete " + file.getAbsolutePath());
	}
	
	private static class Holder
	{
		List<Double> actual;
		List<Double> predicted;
	}
	
	private static Holder runATrial(MetaboliteClass metabolite, String phyla, OtuWrapper wrapper, List<String> keys, 
			boolean scramble) 
			throws Exception
	{
		int otuIndex = wrapper.getIndexForOtuName(phyla);
		HashMap<Integer, List<Double>> metaboliteMap = getMetabolites(metabolite, scramble);
		
		int halfPoint = keys.size() / 2;
		
		File trainingSetFile = new File(ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
				"trainingSet.txt");

		deleteOrThrow(trainingSetFile);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(
								trainingSetFile));
		
		for( int x=0; x < halfPoint; x++)
		{
			int sampleIndex = wrapper.getIndexForSampleName(keys.get(x));
			writer.write( wrapper.getDataPointsNormalizedThenLogged().get(sampleIndex).get(otuIndex)+ " "  );
			
			List<Double> list = metaboliteMap.get( Integer.parseInt(keys.get(x).replace("sample", "")));
			
			for(  int y=0; y < list.size(); y++ )
			{
				boolean skip = false;
				
				if(  metabolite.equals(MetaboliteClass.METADATA) && Math.abs( list.get(y) + 1 ) < 0.0001)
					skip =true;
					
				if( ! skip)
					writer.write( (y+1) + ":" + list.get(y) + " " );
				else
					System.out.println("SKIP!!!!!!!!");
				
			}
				
			writer.write("\n");
		}
		
		writer.flush();  writer.close();
		
		File regressModel = new File(ConfigReader.getMicrboesVsMetabolitesDir() + File.separator +  
				"regressionModel");
		
		deleteOrThrow(regressModel);
		
		// train with svm_learn.exe -z r trainingSet.txt regressModel
		String[] args = new String[5];
		args[0] = ConfigReader.getSvmDir() + File.separator + "svm_learn.exe";
		args[1] = "-z";
		args[2] = "r";
		args[3] = trainingSetFile.getAbsolutePath();
		args[4] = regressModel.getAbsolutePath();
		new ProcessWrapper(args);
		
		File setToClassify = new File( 
				ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
				"setToClassify.txt"	);
		deleteOrThrow(setToClassify);
		
		writer = new BufferedWriter(new FileWriter(setToClassify));
		
		for( int x=halfPoint; x < keys.size(); x++ ) 
		{
			writer.write(keys.get(x) + " " );
			
			List<Double> list = metaboliteMap.get(Integer.parseInt(keys.get(x).replace("sample", "")));
			
			for(  int y=0; y < list.size(); y++ )
				writer.write( (y+1) + ":" + list.get(y) + " " );
			
			writer.write("\n");
		}
		
		writer.flush(); writer.close();
		
		File svmOut = new File(ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
				"svmOut.txt");
		
		deleteOrThrow(svmOut);
		
		// classify with svm_classify setToClassify.txt regressModel svmOut.txt
		args = new String[4];
		args[0] = ConfigReader.getSvmDir() + File.separator + "svm_classify.exe";
		args[1] = setToClassify.getAbsolutePath();
		args[2] = regressModel.getAbsolutePath();
		args[3] = svmOut.getAbsolutePath();
		new ProcessWrapper(args);
		
		BufferedReader svmOutReader = new BufferedReader(new FileReader(svmOut));
		BufferedReader classifiedReader = new BufferedReader(new FileReader(setToClassify));
		List<Double> predicted = new ArrayList<Double>();
		List<Double> actual = new ArrayList<Double>();
		
		for(String s= classifiedReader.readLine(); s != null; s = classifiedReader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			String  sample = sToken.nextToken();
			int sampleID = wrapper.getIndexForSampleName(sample);
			
			actual.add(wrapper.getDataPointsNormalizedThenLogged().get(sampleID).get(otuIndex));
			
			predicted.add(Double.parseDouble(svmOutReader.readLine()));
		}
		
		if( svmOutReader.readLine() != null)
			throw new Exception("No");

		classifiedReader.close();
		svmOutReader.close();
		
		Holder h= new Holder();
		h.actual= actual;
		h.predicted = predicted;
		return h;
		
	}
	
	/*
	 * Scramble will not preserve IDs across categories..
	 */
	private static HashMap<Integer, List<Double>> getMetabolites(MetaboliteClass metaboliteClass,
			boolean scramble) throws Exception
	{

		if( metaboliteClass.equals(MetaboliteClass.BOTH))
		{
			HashMap<Integer, List<Double>> map1 = 
					getMetabolites(MetaboliteClass.PLASMA,scramble);
			
			HashMap<Integer, List<Double>> map2 = 
					getMetabolites(MetaboliteClass.URINE,scramble);
			
			if( ! map1.keySet().equals(map2.keySet()) )
				throw new Exception("Unexpected parse");
			
			for(Integer i : map1.keySet())
			{
				map1.get(i).addAll(map2.get(i));
			}
			
			return map1;
		}
		
		if( metaboliteClass.equals(MetaboliteClass.ALL))
		{
			HashMap<Integer, List<Double>> map1 = 
					getMetabolites(MetaboliteClass.BOTH, scramble);
			
			HashMap<Integer, List<Double>> map2 = 
					getMetabolites(MetaboliteClass.METADATA, scramble);
			
			for(Integer i : map2.keySet())
				if( map1.containsKey(i))
				{
					map1.get(i).addAll(map2.get(i));
				}
				
			return map1;
		}
		
		HashMap<Integer, List<Double>> map = new LinkedHashMap<Integer, List<Double>>();
		
		BufferedReader reader = null;
		
		if( metaboliteClass.equals(MetaboliteClass.URINE))
		{
			reader = 
					new BufferedReader(new FileReader(new File( 
							ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
							"urine_metabolites_data.txt")));
					
		}
		else if ( metaboliteClass.equals(MetaboliteClass.PLASMA))
		{
			reader = 
					new BufferedReader(new FileReader(new File( 
							ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
							"plasmaMetabolites.txt")));
		}
		else if ( metaboliteClass.equals(MetaboliteClass.METADATA))
		{
			reader = 
					new BufferedReader(new FileReader(new File( 
							ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
							"patientMetadataSubjectsAsColumns.txt")));
		}
		else throw new Exception("No");
		
		String firstLine = reader.readLine();
		
		TabReader tr = new TabReader(firstLine);
		
		tr.nextToken();
		
		if (! metaboliteClass.equals(MetaboliteClass.METADATA))
		{
			tr.nextToken(); tr.nextToken(); 
		}
		
		while( tr.hasMore() )
		{
			String nextToken = tr.nextToken();
			
			if( nextToken.trim().length() > 0 )
			{
				int aVal = Integer.parseInt(nextToken);
				map.put(aVal, new ArrayList<Double>());
			}
		}
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			tr =new TabReader(s);
			
			tr.nextToken(); 
		
			if (! metaboliteClass.equals(MetaboliteClass.METADATA))
			{
				tr.nextToken(); tr.nextToken(); 
			}
			
			for( Integer key : map.keySet())
			{
				List<Double> list =map.get(key);
				list.add(Double.parseDouble(tr.nextToken()));
			}

			if( tr.hasMore() && tr.nextToken().trim().length() > 0 )
				throw new Exception("No " + tr.nextToken());
		}
		
		reader.close();
		
		if( scramble)
		{
			HashMap<Integer, List<Double>> newMap = new LinkedHashMap<Integer, List<Double>>();
			List<Integer> newKeys = new ArrayList<Integer>(map.keySet());
			Collections.shuffle(newKeys);
			
			int counter = 0;
			
			for( List<Double> innerList : map.values() )
			{
				newMap.put(newKeys.get(counter), innerList);
				counter++;
			}
			
			if( map.size() != newMap.size())
				throw new Exception("No");

			if( map.keySet().size() != newMap.keySet().size())
				throw new Exception("No");
			
			if( ! map.keySet().equals( newMap.keySet()))
				throw new Exception("No");
			
			map = newMap;
				
		}
		
		for( List<Double> list : map.values() )
		{
			double sum =0;
			int n=0;
			
			for( int y=0; y < list.size(); y++)
			{
				boolean addin = true;
				
				if( metaboliteClass.equals(MetaboliteClass.METADATA) && 
					Math.abs(	list.get(y) + 1 ) <= 0.0001 ) 
					addin = false;
				
				if( addin)
				{
					sum += list.get(y);
					n++;
				}
			}
			
			if ( n > 0 )
			{
				double average = sum / n;
				for( int y=0; y < list.size(); y++)
				{
					list.set(y, list.get(y) / average);
				}
			}			
		}
		
		return map;
	}
	
	public static void main(String[] args) throws Exception
	{
		OtuWrapper wrapper = new OtuWrapper(ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + "phylaAsColumns.txt");
		
		for( String s : wrapper.getOtuNames())
		{
			writeATrialFile(wrapper, s, false);
			writeATrialFile(wrapper, s, true);
		}
	}
	 
	public static void writeATrialFile(OtuWrapper wrapper, String phyla, boolean scramble) throws Exception
	{
		Random random= new Random(324234);
		
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File( 
			ConfigReader.getMicrboesVsMetabolitesDir() + File.separator + 
			"trials_comp"+ phyla + "_scrambled_" +scramble +  ".txt")));
		
		writer.write("pValuePlasma\trValuePlasma\tpValueUrine\trValueUrine\tpValueBoth\trValueBoth\t");
		writer.write("pValueMetadata\trValueMetadata\tpValueAll\trValueAll\n");
		
		MetaboliteClass mClass[] = { MetaboliteClass.PLASMA, MetaboliteClass.URINE, MetaboliteClass.BOTH,
				MetaboliteClass.METADATA, MetaboliteClass.ALL};
		
		for( int x=0; x < 10; x++)
		{
			List<String> keys = new ArrayList<String>(wrapper.getSampleNames());
			Collections.shuffle(keys, random);
			for( int y=0; y < mClass.length; y++)
			{
				Holder h = runATrial(mClass[y], phyla,wrapper, keys, scramble);
				Regression r = new Regression();
				r.fitFromList(h.actual, h.predicted);
				writer.write(r.getPValueForSlope()+ "\t");
				writer.write( Pearson.getPearsonR(h.actual, h.predicted) + 
						(y==4 ? "\n" : "\t") );
				writer.flush();
			}
		}
		
		writer.flush();  writer.close();
	}
}
