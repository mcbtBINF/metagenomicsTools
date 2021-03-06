package scripts.goranLab;

import java.io.File;

import parsers.OtuWrapper;
import utils.ConfigReader;

public class RePivotOTU
{
	public static void main(String[] args) throws Exception
	{
		OtuWrapper.transpose(ConfigReader.getGoranTrialDir() + File.separator + 
				"otuCounts.txt",ConfigReader.getGoranTrialDir() + File.separator + 
				File.separator + "otu" + "fromOTUsAsColumn.txt" , false);
		
		OtuWrapper wrapper = new OtuWrapper(ConfigReader.getGoranTrialDir() + File.separator + 
				File.separator + "otu" + "fromOTUsAsColumn.txt");
		
		wrapper.writeLoggedDataWithTaxaAsColumns(new File(ConfigReader.getGoranTrialDir() + File.separator + 
				File.separator + "otu" + "fromOTUsAsColumnLogNorm.txt"));
	}
}
