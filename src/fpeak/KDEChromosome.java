package fpeak;
import java.util.Arrays;
import java.util.Random;
import java.io.File;

import background.BffReader;
import background.IffReader;
import background.WigChromosome;
import io.DensityWriter;


public class KDEChromosome {
  
  static int BATCH_SIZE = 1024 * 10;
  
  private long _firstCut;
  private long _lastCut;
  private Sequence[] _cuts;
  private String _chromosome;
  private float _threshold;
  private int _sequenceLength;
  
  public KDEChromosome(String chromosome, Sequence[] cuts, int sequenceLength){
    _chromosome = chromosome;
    _cuts = cuts;
    _firstCut = cuts[0].getPosition();
    _lastCut = cuts[cuts.length-1].getPosition();
    _sequenceLength = sequenceLength;
  }
  
  public int getSequenceLength() {
	  return _sequenceLength;
  }
  
  public Sequence[] getCuts() {
	  	return _cuts;
  }
  
  public String getChromosome(){
    return _chromosome;
  }
  
  public long getFirstPos(){
    return _firstCut;
  }
  
  public long getLastPos(){
    return _lastCut;
  }
  
  public long getLength(){
	    return _cuts.length;
  }
  
  /*
   * This function is the same as above but reads in a background uniqueness file and adjusts the
   * density estimate based on this.
   */
  public void runBG(Settings settings, DensityWriter dw, boolean verboseFlag, float wg_threshold, File bgfile) throws Exception {
//    _threshold = computeThreshold(settings);
	    _threshold = wg_threshold;
	    dw.setThreshold(_threshold);
	    
		WigChromosome bgchr = null;
		try {
			bgchr = BffReader.read(bgfile);
		} catch (Exception e){
			e.printStackTrace();
		}
//		System.out.println(chr.getChromosome());
//		System.out.println(chr.getLength());
//		System.out.println(chr.getStart());
//	    if(verboseFlag){
//	    	System.out.println(_chromosome + ": Threshold: " + _threshold);
//	    }
	    
	    Sequence[] cuts = _cuts;
//	    long[] cuts = _cuts;
	    
	    float[] density = new float[BATCH_SIZE];
	    
	    if(verboseFlag){
	      System.out.println(_chromosome + ": first=" + _firstCut + ", last=" + _lastCut);
	      for(int i = 0; i < 20; ++i){
	        System.out.print(".");
	      }
	      System.out.println();
	    }
	    int numBases = (int)Math.abs(_lastCut - _firstCut);
	    
	    int incr =  numBases / 20;
	    int mod = BATCH_SIZE -1;
	    int peaks = 0;
	    
	    long start = System.currentTimeMillis();
	    int cutIdx = 0;
	    long currentChromPos = 0;
	    int arrPos = 0;
	    boolean aboveThreshold = false;
	    for(int i = 0; i < numBases; ++i){
	      currentChromPos = i + _firstCut;
	      arrPos = i % BATCH_SIZE;
	      density[arrPos] = (float)bgdensity(settings, currentChromPos, cutIdx, cuts, bgchr);
	      
	      if(!aboveThreshold && density[arrPos] > _threshold){
	        aboveThreshold = true;
	        ++peaks;
	      }else if(aboveThreshold && density[arrPos] < _threshold){
	        aboveThreshold = false;
	      }
	      while(cutIdx < cuts.length && currentChromPos > cuts[cutIdx].getPosition())
	        ++cutIdx;
	      
	      if(verboseFlag && i % incr == 0)
	        System.out.print(".");
	      if(i % BATCH_SIZE == mod){
	        dw.writeDensity(density, 0, BATCH_SIZE); 
	      }
	    }
	    
	    int len = numBases % BATCH_SIZE;
	    dw.writeDensity(density, 0, len);
	    
	    if(verboseFlag){
	      System.out.println();
	      System.out.println(_chromosome + ": Completed in " + (System.currentTimeMillis() - start)/1000d + " seconds.");
	      System.out.println(_chromosome + ": Found " + peaks + " peaks.");
	      System.out.println("-----------------------");
	    }
  }

  /*
   * This function is the same as above but reads in a background uniqueness file and adjusts the
   * density estimate based on this.
   */
  public void run(Settings settings, DensityWriter dw, boolean verboseFlag, float wg_threshold, File[] bgfile, File[] ipfile) throws Exception {
//    _threshold = computeThreshold(settings);
	    _threshold = wg_threshold;
	    dw.setThreshold(_threshold);
//	    dw.setThreshold(10);
	    
	    boolean bg_hit = false;
	    boolean ip_hit = false;
	    boolean bg_used = true;
	    boolean ip_used = true;
	    File bg_file_used = null;
	    File ip_file_used = null;
	    
	    if(bgfile.length == 0) {
	    	// No background used
	    	bg_used = false;
	    } else {
	    	for(int j = 0; j < bgfile.length; ++j) {
	    		if(bgfile[j].getName().equals(_chromosome + ".bff")) {
	    			bg_file_used = bgfile[j];
	    			//System.out.println("Running background on Chromosome " + _chromosome);
	    			bg_hit = true;
	    		}
	    	}
	    	if(!bg_hit) {
	    		System.out.println("No background for Chromosome " + _chromosome);
	    		return;
	    	}
	    }
	    
	    if(ipfile.length == 0) {
	    	// No ploidy used
	    	ip_used = false;
	    } else {
	    	for(int j = 0; j < ipfile.length; ++j) {
	    		if(ipfile[j].getName().equals(_chromosome + ".iff")) {
	    			ip_file_used = ipfile[j];
	    			//System.out.println("Running input on Chromosome " + _chromosome);
	    			ip_hit = true;
	    		}
	    	}
	    	if(!ip_hit) {
	    		System.out.println("No input for Chromosome " + _chromosome);
	    		return;
	    	}
	    }
	    
	    //Load files
    	WigChromosome bgchr = null;
    	WigChromosome ipchr = null;
	    if(bg_used) {
	    	try {
	    		bgchr = BffReader.read(bg_file_used);
	    	} catch (Exception e){
	    		e.printStackTrace();
	    	}
	    }
	    if(ip_used) {
	    	try {
	    		ipchr = IffReader.read(ip_file_used);
	    	} catch (Exception e){
	    		e.printStackTrace();
	    	}
	    }
//		System.out.println(chr.getChromosome());
//		System.out.println(chr.getLength());
//		System.out.println(chr.getStart());
//	    if(verboseFlag){
//	    	System.out.println(_chromosome + ": Threshold: " + _threshold);
//	    }
	    
	    Sequence[] cuts = _cuts;
//	    long[] cuts = _cuts;
	    
	    float[] density = new float[BATCH_SIZE];
	    float[] densityP = new float[BATCH_SIZE];
	    float[] densityM = new float[BATCH_SIZE];
	    
	    if(verboseFlag){
	      System.out.println(_chromosome + ": first=" + _firstCut + ", last=" + _lastCut);
	      for(int i = 0; i < 20; ++i){
	        System.out.print(".");
	      }
	      System.out.println();
	      System.out.println("running the P/M version");
	    }
	    int numBases = (int)Math.abs(_lastCut - _firstCut);
	    
	    int incr =  numBases / 20;
	    int mod = BATCH_SIZE -1;
	    int peaks = 0;
	    
	    long start = System.currentTimeMillis();
	    int cutIdx = 0;
	    long currentChromPos = 0;
	    int arrPos = 0;
	    boolean aboveThreshold = false;
	    float[] sumPM=new float[2];
	    for(int i = 0; i < numBases; ++i){
	      currentChromPos = i + _firstCut;
	      arrPos = i % BATCH_SIZE;
	      
	      if(!bg_used && !ip_used) {
	    	  density[arrPos] = (float)densityPM(settings, currentChromPos, cutIdx, cuts, sumPM);
	    	  densityP[arrPos] = sumPM[0];
	    	  densityM[arrPos] = sumPM[1];
	      } else {
	    	  if(bg_used && !ip_used) {
//	    		  density[arrPos] = (float)bgdensity(settings, currentChromPos, cutIdx, cuts, bgchr);
	    		  density[arrPos] = (float)bgdensityPM(settings, currentChromPos, cutIdx, cuts, bgchr, sumPM);
	    		  densityP[arrPos] = sumPM[0];
		    	  densityM[arrPos] = sumPM[1];
	    	  } else {
	    		  if(!bg_used && ip_used) {
	    			  density[arrPos] = (float)ipdensity(settings, currentChromPos, cutIdx, cuts, ipchr);
	    		  } else {
//	    			  density[arrPos] = (float)fulldensity(settings, currentChromPos, cutIdx, cuts, bgchr, ipchr);
	    			  density[arrPos] = (float)fulldensityPM(settings, currentChromPos, cutIdx, cuts, bgchr, ipchr, sumPM);
	    			  densityP[arrPos] = sumPM[0];
			    	  densityM[arrPos] = sumPM[1];
	    		  }
	    	  }
	      }
//	    	  density[arrPos] = (float)bgdensity(settings, currentChromPos, cutIdx, cuts, bgchr);
	      
	      if(!aboveThreshold && density[arrPos] > _threshold){
	        aboveThreshold = true;
	        ++peaks;
	      }else if(aboveThreshold && density[arrPos] < _threshold){
	        aboveThreshold = false;
	      }
	      while(cutIdx < cuts.length && currentChromPos > cuts[cutIdx].getPosition())
	        ++cutIdx;
	      
	      if(verboseFlag && i % incr == 0)
	        System.out.print(".");
	      if(i % BATCH_SIZE == mod){
//	        dw.writeDensity(density, 0, BATCH_SIZE);
	    	  dw.writeDensityPM(density, 0, BATCH_SIZE, densityP, densityM);
	      }
	    }
	    
	    int len = numBases % BATCH_SIZE;
//	    dw.writeDensity(density, 0, len);
	    dw.writeDensityPM(density, 0, len, densityP, densityM);
	    
	    if(verboseFlag){
	      System.out.println();
	      System.out.println(_chromosome + ": Completed in " + (System.currentTimeMillis() - start)/1000d + " seconds.");
	      System.out.println(_chromosome + ": Found " + peaks + " peaks.");
	      System.out.println("-----------------------");
	    }
  }
  /*  
  private float computeThreshold(Settings settings){
    Random r = new Random();
    
    double size = (int)Math.abs(_lastCut - _firstCut);
    double ncuts = _cuts.length;
    int totalWindow = 1 + (int)(settings.window * 2);
    
    int cutDensity = (int)((ncuts / size) * totalWindow);
    int thresholdIterations = 1000;
    long[] cuts = new long[cutDensity];
    double[] densities = new double[thresholdIterations];
    
    
    for(int i = 0; i < thresholdIterations; ++i){
      for(int j = 0; j < cuts.length; ++j)
        cuts[j] = r.nextInt(totalWindow);
      Arrays.sort(cuts);
      //densities[i] = density(settings, (long)settings.window, (int)(cuts.length/2.0), cuts);
    }
    
    double mean = Util.mean(densities);
    double std = Util.std(densities);
    
    return (float)(mean + settings.threshold * std);
  }
*/
  
  private static float density(Settings settings, long chromPos, int cutIdx, Sequence[] cuts){
    
    long minPos = chromPos - settings.window;
    long maxPos = chromPos + settings.window;
    
    double[] PRECOMPUTE = settings.precompute;
    
    double sum = 0.0;
    double sumP=0.0, sumM=0.0;
    for(int i = cutIdx-1; i > -1; --i){
      if (cuts[i].getPosition() < minPos) 
        break;
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      if(!settings.dnaseExperimentType) {
    	  if (cuts[i].getStrand())
    		  sumP+=settings.precompute[d];
    	  else sumM+=settings.precompute[d];
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  d = Math.abs((int)(cuts[i].getPosition() + settings.offset - chromPos));
    		  sum += settings.precompute[d];
    	  } else {
    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
    			  d = Math.abs((int)(cuts[i].getPosition() - settings.offset - chromPos));
    			  sum += settings.precompute[d];
    		  }
    	  }
      } else {
    	  sum += settings.precompute[d];
      }
    }
    
    for(int i = cutIdx; i < cuts.length; ++i){
      if (cuts[i].getPosition() > maxPos) break;
      
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      
      //System.out.println(d);
      if(d > PRECOMPUTE.length-1)
        throw new IllegalStateException();
      
      if(!settings.dnaseExperimentType) {
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  sumP+=settings.precompute[d];
    		  d = Math.abs((int)(cuts[i].getPosition() + settings.offset - chromPos));
    		  sum += settings.precompute[d];
    	  } else {
    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
    			  sumM+=settings.precompute[d];
    			  d = Math.abs((int)(cuts[i].getPosition() - settings.offset - chromPos));
    			  sum += settings.precompute[d];    		  
    		  }
    	  }
      } else {
    	  sum += settings.precompute[d];
      }
    }
    
    return (float)(sum / (double)settings.bandwidth);
//    return (float)Math.log(sumP/sumM);
  }
  
  private static float densityPM(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, float[] sumPM){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	    
	    double sum = 0.0;
	    sumPM[0]=0.0f;
	    sumPM[1]=0.0f;
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      if(!settings.dnaseExperimentType) {
	    	  if (d<20){
		    	  if (cuts[i].getStrand())
		    		  sumPM[0]+=settings.precompute[d];
		    	  else sumPM[1]+=settings.precompute[d];
	    	  }
	    	  if(cuts[i].getStrand()) {
//	    		  d = Math.abs((int)(cuts[i].getPosition() + settings.offset - chromPos));
	    		  sum += settings.precompute[d];
	    	  } else {
//	    			  d = Math.abs((int)(cuts[i].getPosition() - settings.offset - chromPos));
	    			  sum += settings.precompute[d];
	    		  }
	      } else {
	    	  sum += settings.precompute[d];
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      //System.out.println(d);
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(!settings.dnaseExperimentType) {
	    	  if (d<20){
		    	  if (cuts[i].getStrand())
		    		  sumPM[0]+=settings.precompute[d];
		    	  else sumPM[1]+=settings.precompute[d];
	    	  }
	    	  if(cuts[i].getStrand()) {
//	    		  d = Math.abs((int)(cuts[i].getPosition() + settings.offset - chromPos));
	    		  sum += settings.precompute[d];
	    	  } else {
//	    			  d = Math.abs((int)(cuts[i].getPosition() - settings.offset - chromPos));
	    			  sum += settings.precompute[d];    		  
	    		  }
	      } else {
	    	  sum += settings.precompute[d];
	      }
	    }
	    
	    return (float)(sum / (double)settings.bandwidth);
//	    return (float)Math.log(sumP/sumM);
	  }
  
  
  private float bgdensity(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata){
    
    long minPos = chromPos - settings.window;
    long maxPos = chromPos + settings.window;
    
    double[] PRECOMPUTE = settings.precompute;
   
    int bgOffset = (int)_firstCut - bgdata.getStart(); //convert to 0 based
    
    double sum = 0.0;
    int b;
    
    for(int i = cutIdx-1; i > -1; --i){
      if (cuts[i].getPosition() < minPos) 
        break;
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      
      
      if(!settings.dnaseExperimentType) {
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
    			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
    		  }
    	  } else {
    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
        			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
        		  }  		  
    		  }
    	  }
      } else {
    	  if(cuts[i].getStrand()) {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
    	  } else {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
    	  }
		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
		  }
    	  //sum += settings.precompute[d];
      }
    }
    
    for(int i = cutIdx; i < cuts.length; ++i){
      if (cuts[i].getPosition() > maxPos) break;
      
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      
      //System.out.println(d);
      if(d > PRECOMPUTE.length-1)
        throw new IllegalStateException();
      
      if(!settings.dnaseExperimentType) {
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
    			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
    		  }
    	  } else {
    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
        			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
        		  }		  
    		  }
    	  }
      } else {
    	  if(cuts[i].getStrand()) {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
    	  } else {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
    	  }
		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
		  }
      }
    }
    
    return (float)(sum / (double)settings.bandwidth);
  }
  
  private float bgdensityPM(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata, float[] sumPM){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	    double sum = 0.0, sig;
	    int b;
	    sumPM[0]=0.0f;
	    sumPM[1]=0.0f;
	    
	    int cntL=0, cntR=0, sumD=0;
	    
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      if (cuts[i].fragLen<200) continue;
	      int maxD=cuts[i].fragLen;
	      if (maxD==0)
	    	  maxD=2*settings.offset;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      if (chromPos==settings.pos){
		    	cntL++;
		    	sumD+=d;
	      }
	      
	      
	      if(!settings.dnaseExperimentType) {
	    	  if (d<20){
		    	  if (cuts[i].getStrand())
		    		  sumPM[0]+=settings.precompute[d];
		    	  else sumPM[1]+=settings.precompute[d];
	    	  }
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
//	    	  if(cuts[i].getStrand()) {
	    		  if (cuts[i].getPosition() <= chromPos){
	    			  if (d<maxD)
	    				  sig=1.0;
	    			  else sig=0;
	    		  }
	    		  else sig=settings.precompute[d];
//	    		  if (maxD>0)
//	    			  sig/=maxD;
//	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)maxD/2 - chromPos));
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0 &&d<settings.window) {
	    			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
//	    			  sum += sig*(double)bgdata.getValues()[b];
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
//    			  if(!cuts[i].getStrand()) {
    				  if (cuts[i].getPosition() >= chromPos){
		    			  if (d<maxD)
		    				  sig=1.0;
		    			  else sig=0;
		    		  }
		    		  else sig=settings.precompute[d];
//    				  if (maxD>0)
//    	    			  sig/=maxD;
//	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
    				  d = Math.abs((int)(cuts[i].getPosition() - (int)maxD/2 - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0 &&d<settings.window) {
	        			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
//	        			  sum += sig*(double)bgdata.getValues()[b];
	        		  }
	    		  }
	    		  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    	  }
			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
				  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
			  }
	    	  //sum += settings.precompute[d];
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      if (cuts[i].fragLen<200) continue;
	      int maxD=cuts[i].fragLen;
	      if (maxD==0)
	    	  maxD=2*settings.offset;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      if (chromPos==settings.pos){
		    	cntR++;
		    	sumD+=d;
	      }
	      
	      //System.out.println(d);
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(!settings.dnaseExperimentType) {
	    	  if (d<20){
		    	  if (cuts[i].getStrand())
		    		  sumPM[0]+=settings.precompute[d];
		    	  else sumPM[1]+=settings.precompute[d];
	    	  }
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
//    		  if(cuts[i].getStrand()) {
    			  if (cuts[i].getPosition() <= chromPos){
	    			  if (d<maxD)
	    				  sig=1.0;
	    			  else sig=0;
	    		  }
	    		  else sig=settings.precompute[d];
//    			  if (maxD>0)
//	    			  sig/=maxD;
//	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
    			  d = Math.abs((int)(cuts[i].getPosition() + (int)maxD/2 - chromPos));
	    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0 &&d<settings.window) {
	    			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
//	    			  sum+=sig* (double)bgdata.getValues()[b];
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
//	    		  if(!cuts[i].getStrand()) {
	    			  if (cuts[i].getPosition() >= chromPos){
		    			  if (d<maxD)
		    				  sig=1.0;
		    			  else sig=0;
		    		  }
		    		  else sig=settings.precompute[d];
//	    			  if (maxD>0)
//		    			  sig/=maxD;
//	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)maxD/2 - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0 &&d<settings.window) {
	        			  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
//	    				  sum+=sig * (double)bgdata.getValues()[b];
	        		  }
	    		  }
	    		  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    	  }
			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
				  sum += settings.precompute[d] * (double)bgdata.getValues()[b];
			  }
	      }
	    }
	    
	    if (chromPos==settings.pos){
	    	// for debugging purpose
	    	System.out.println(chromPos+"\t"+cntL+"\t"+cntR+"\t"+sumD/(cntL+cntR)+"\t"+sum / (double)settings.bandwidth);
	    	System.exit(0);
	    }
	    return (float)(sum / (double)settings.bandwidth);
//	    return (float)sum;
	  }
  
  private float ipdensity(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	   
	    int bgOffset = (int)_firstCut - bgdata.getStart(); //convert to 0 based
	    
	    double sum = 0.0;
	    int b;
	    
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      
	      if(!settings.dnaseExperimentType) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d);
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	        			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d);
	        		  }  		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    	  }
			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
				  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d);
			  }
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      //System.out.println(d);
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(!settings.dnaseExperimentType) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d);
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	        			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d);
	        		  }		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    	  }
			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
				  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d);
			  }
	      }
	    }
	    
	    return (float)(sum / (double)settings.bandwidth);
	  }
  
  private float fulldensity(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata, WigChromosome ipdata){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	   
	    int bgOffset = (int)_firstCut - bgdata.getStart(); //convert to 0 based
	    
	    double sum = 0.0;
	    int b;
	    int c;
	    
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      
	      if(!settings.dnaseExperimentType) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart(); //index of ip for particular sequence i
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    			  }
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
		    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
		    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
		    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
		    			  }
	        		  }  		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    	  }
    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
    			  }
    		  }  		
	    	  //sum += settings.precompute[d];
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      //System.out.println(d);
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(!settings.dnaseExperimentType) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000);
	    			  }	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
		    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
		    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
		    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
		    			  }	        		  }		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    	  }
    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
    			  }
    		  }  		
	    	  //sum += settings.precompute[d];
	      }
	    }
	    
	    return (float)(sum / (double)settings.bandwidth);
	  }
  
  private float fulldensityPM(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata, WigChromosome ipdata, float[] sumPM){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	    
	    double sum = 0.0, sig;
	    int b;
	    int c;
	    sumPM[0]=0.0f;
	    sumPM[1]=0.0f;
	    
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      
	      if(!settings.dnaseExperimentType) {
	    	  if (d<20){
		    	  if (cuts[i].getStrand())
		    		  sumPM[0]+=settings.precompute[d];
		    	  else sumPM[1]+=settings.precompute[d];
	    	  }
	    	  if(cuts[i].getStrand()) {
	    		  if (cuts[i].getPosition() <= chromPos){
	    			  if (d<2*settings.offset)
	    				  sig=1.0;
	    			  else sig=0;
	    		  }
	    		  else sig=settings.precompute[d];
//	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart(); //index of ip for particular sequence i
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
//	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    				  sum += sig * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    			  }
	    		  }
	    	  } else {
	    		  if (cuts[i].getPosition() >= chromPos){
	    			  if (d<2*settings.offset)
	    				  sig=1.0;
	    			  else sig=0;
	    		  }
	    		  else sig=settings.precompute[d];
//	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
//	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    				  sum += sig * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    			  }
        		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    	  }
  		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
  			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
  				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
  			  }
  		  }  		
	    	  //sum += settings.precompute[d];
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      //System.out.println(d);
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(!settings.dnaseExperimentType) {
	    	  if (d<20){
		    	  if (cuts[i].getStrand())
		    		  sumPM[0]+=settings.precompute[d];
		    	  else sumPM[1]+=settings.precompute[d];
	    	  }
	    	  if(cuts[i].getStrand()) {
	    		  if (cuts[i].getPosition() <= chromPos){
	    			  if (d<2*settings.offset)
	    				  sig=1.0;
	    			  else sig=0;
	    		  }
	    		  else sig=settings.precompute[d];
//	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
//	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000);
	    				  sum += sig * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000);
	    			  }	    		  }
	    	  } else {
	    		  if (cuts[i].getPosition() >= chromPos){
	    			  if (d<2*settings.offset)
	    				  sig=1.0;
	    			  else sig=0;
	    		  }
	    		  else sig=settings.precompute[d];
//	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
//		    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    				  sum += sig * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
	    			  }
    			  }
    		  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    	  }
  		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
  			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
  				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d);
  			  }
  		  }  		
	    	  //sum += settings.precompute[d];
	      }
	    }
	    
//	    return (float)(sum / (double)settings.bandwidth);
	    return (float)sum;
	  }
  
  public static class Sequence implements Comparable
  {
      long position;
      boolean strand;
      int fragLen;
      
      public long getPosition() {
    	  	return position;
      }
      
      public boolean getStrand() {
    	  	return strand;
      }
      
      public void setPosition(long pos) {
    	  this.position = pos;
      }
    
      public void setStrand(boolean str) {
    	  this.strand = str;
      }
      
      public Sequence(long x, boolean y, int len) {
          this.position = x;
          this.strand = y;
          this.fragLen= len;
      }
      
      public Sequence() {
    	  this.position = 0;
    	  this.strand = false;
    	  this.fragLen= 0;
      }
      
      //CompareTo overload
      public int compareTo(Object obj) {
    	  Sequence tmp = (Sequence)obj;
    	  if(this.position < tmp.position) {
    		  return -1;
    	  }else if(this.position > tmp.position) {
    		  return 1;
    	  }
    	  return 0;
      }
  }

  public static class Settings {
    private static final double PI2 = Math.sqrt(Math.PI * 2);
    static final float DEFAULT_THRESHOLD = 4.0f;
    
    public final double[] precompute;
    public final long window;
    public final long bandwidth;
    public final int step;
    public final float threshold;
    public final int offset;
    public final boolean dnaseExperimentType;
    public final long ncuts;
    public long pos;
    
    public Settings(long bandwidth, long window, float threshold, int offset, long ncuts){
      this.bandwidth = bandwidth;
      this.window = window;
      this.threshold = threshold;
      this.offset = offset;
      this.step = 0;
      this.ncuts = ncuts;
      
      dnaseExperimentType = isDNase(offset);
      precompute = precompute(window, bandwidth, ncuts);
    }
    
    public Settings(long featureLength, float threshold, int offset, long ncuts, long pos){
      this.bandwidth = computeBandwidth(featureLength);
      this.window = computeOptimalWindow(bandwidth);
      this.threshold = threshold;
      this.offset = offset;
      this.step = 0;
      this.ncuts = ncuts;
      this.pos=pos;
      
      dnaseExperimentType = isDNase(offset);
      precompute = precompute(window, bandwidth, ncuts);
    }
    
    private static double sequenceNormalize(double value, long ncuts) {
  	  return ((value * 20000000d)/(double)ncuts);
    }//Add this to all density outputs.
    
    private static double[] precompute(long window, long bandwidth, long ncuts){
      double[] precompute = new double[(int)(window + 1)];
      for(int i = 0; i < precompute.length; ++i){
        double x = i / (double)bandwidth;
        precompute[i] = sequenceNormalize( Math.exp(-(x * x) / 2) / PI2 , ncuts);
//        precompute[i] =  Math.exp(-(x * x) / 2);
      }
      return precompute;
    }
    
    private static boolean isDNase(int offset){
    	if(offset == 0) {
    		return true;
    	}
    	return false;
    }
    
    public static long computeBandwidth(long featureLength){
      return (long)((featureLength / 2.0) / 3.0);  // 3 standard deviations
    }

    public static int computeOptimalWindow(long bandwidth){
      int i = 1400;
      
      double bw = (double)bandwidth;
      while(true){
        double x = ++i / bw;
        double v = Math.exp(-(x * x) / 2) / PI2;
        if(v < Float.MIN_VALUE)
          break;
      }
      return (i-1);
    }
  }
  
}
