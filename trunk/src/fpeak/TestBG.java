package fpeak;
import java.io.File;
import java.io.FileFilter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import background.IffReader;
import background.WigChromosome;

public class TestBG {
	
	public static void main(String args[]){
		WigChromosome ipchr = null;
		File ip_file_used = null;
		Options opts = new Options();
		String ploidyDirectory = null;
		String[] ipfiles = {};
		File[] ploidy_files=null;
		opts.addOption(OptionBuilder.withArgName( "ploidy dir" )
				.hasArg()
	            .withDescription( "ploidy/input directory (default=none)" )
	            .isRequired(false)
	            .create( "p" ));
		CommandLineParser parser = new GnuParser();
		try{
			CommandLine cmd = parser.parse(opts, args);
			if(cmd.hasOption("p")) //ploidy|input directory
		        ploidyDirectory = cmd.getOptionValue("p");
			File[] ploidy = getFiles(ploidyDirectory, ipfiles);
			ipchr = IffReader.read(ploidy[0]);
			for (int i=0;i<ipchr.getLength();i++){
				System.out.println((i+ipchr.getStart())+": "+ipchr.getValues()[i]);
			}
		}
		catch (Exception ex){
			ex.printStackTrace();
		}
	}
	
	private static File[] getFiles(String directory, String[] files){
		  
	    File[] inputFiles = null;
	    
	    if(directory != null){
	      File dir = new File(directory);
	      if (!dir.exists() || !dir.isDirectory()){
	        System.out.println("Directory parameter " + directory + " does not exist or is not a directory.");
	        System.exit(1);
	      }
	      if(files.length == 0){
	        inputFiles = dir.listFiles(new FileFilter(){
	          public boolean accept(File f) {
	            return !f.isHidden() && f.isFile();
	          }
	        });
	      }
	    }else{
	      directory = System.getProperty("user.dir");
	    }
	    
	    if(inputFiles == null){
	      inputFiles = new File[files.length];
	      for(int i = 0; i < files.length; ++i){
	        File f = new File(files[i]);
	        boolean exists = false;
	        if(!(exists = f.exists())){
	          f = new File(directory, files[i]);
	          exists = f.exists();
	        }
	        if(!exists){
	          System.out.println("Input file " + files[i] + " does not exist.");
	          System.exit(1);
	        }
	        inputFiles[i] = f;
	      }
	    }
	   
	    return inputFiles;
	  }
}