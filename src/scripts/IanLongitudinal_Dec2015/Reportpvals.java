import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

public class Reportpvals {
    public static void main(String[] args){
	String dirName = "/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal/";
	String[] taxaLevels = {"phylum", "class", "order", "family", "genus"};
	// Grab each of the model files.
	// Should probably read list of files I care about from a file.
	String[] modelBaseFileName = {"pValuesLongPatient_BMI_ANOVA_NoLow","pValuesLongPatient_Day_ANOVA_NoLow"," pValuesLongPatient_EnergyIntake_ANOVA_NoLow"};
	File[] pVals = finder(dirName, modelBaseFileName);
	for(int i = 0; i < pVals.length; i++){
	    //Open up a file to write to
	    for(int j = 0; j < 5; j++){
		String specificName = pvals[i].toString() + "_" + taxaLevels[j];
		writeContents(pvals[i], taxaLevels[j], readFile((File) specificName));
	    }
	    //Close that file
	}
    }

    public static File[] finder(String dirName, String[] modelNames){
	File dir = new File(dirName);

	return dir.listFiles(new FilenameFilter()
	{
	    public boolean accept(File dir, String filename)
    	    {
		//Go over the model names too.
		return filename.endsWith("fastatoRDP.txt");
	    }
	} );
    }

    public static void writeContents(File fout, String preString, String[] preprocessed){

    }

    //This is probably a standard thing done in java.
    public static String[] readFile(File f){
	//Open the file
	String[] fileasStringArray;
	BufferedReader reader = new BufferedReader(new FileReader(f));
	String nextLine = reader.readLine();
	while(nextLine != null){
	    fileasStringArray
	}

	return fileasStringArray;
    }
}
