package io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import java.text.NumberFormat;

import fpeak.KDEChromosome;
import fpeak.KDEChromosome.Settings;
import background.BffReader;
import background.IffReader;
import background.WigChromosome;

public class NpfDensityWriter implements DensityWriter{

  private BufferedWriter bw;
  private final long chromPos;
  private final int step;
  private final String chr;
  
  private float _threshold = 0.0f;
  private long _currentPos;
  private long _counter = 1;
  
  private long _startPeakPos = 0;
  private boolean _aboveThreshold = false;
  private float _currentMax = 0.0f;
  private float _currentMean = 0.0f;
  private long _currentMaxPos = 0;
  
  private NumberFormat nf;
  
  public NpfDensityWriter(File f, String chr, long chromStart, int step) throws IOException {
	    bw = new BufferedWriter(new FileWriter(f));
	    this.chromPos = chromStart;
	    this.chr = chr;
	    this.step = step;
	    nf = NumberFormat.getNumberInstance();
	    nf.setGroupingUsed(false);
	    nf.setMaximumFractionDigits(8);
	    nf.setMinimumFractionDigits(8);
	    
	    _currentPos = chromPos;
	  }
	  
	  public void setThreshold(float threshold){
	    _threshold = threshold;
	  }
	  
	  public void close() throws IOException {
	    if(_aboveThreshold)
	      doWrite();
	    bw.close();
	  }

	  public void writeDensity(float[] batch, int start, int length)
	      throws IOException {
	    int end = start + length;
	    for(int i = start; i < end; ++i, ++_currentPos){
	      if(!_aboveThreshold){
	        if(batch[i] > _threshold){
	          _aboveThreshold = true;
	          _startPeakPos = _currentPos;
	          _currentMax = batch[i];
	          _currentMaxPos = _currentPos;
	        }
	      }else{ // aboveThreshold
	        if(batch[i] > _threshold){
	          _currentMax = Math.max(_currentMax, batch[i]);
	          if(_currentMax == batch[i])
	        	  _currentMaxPos = _currentPos;
	        }else{
	          _aboveThreshold = false;
	          doWrite();
	        }
	      }
	    }
	  }
	  
	  public void writeDensityPM(float[] batch, int start, int length, float[] batchP, float[] batchM)
		      throws IOException {
		    int end = start + length;
		    int cent=0,beg=0;
		    for(int i = start; i < end; ++i, ++_currentPos){
		      if(!_aboveThreshold){
		        if(batch[i] > _threshold){
		          _aboveThreshold = true;
		          _startPeakPos = _currentPos;
		          _currentMax = batch[i];
		          _currentMaxPos = _currentPos;
		          _currentMean = batch[i];
		          cent=i;
		          beg=i;
		        }
		      }else{ // aboveThreshold
		        if(batch[i] > _threshold){
		          _currentMax = Math.max(_currentMax, batch[i]);
		          _currentMean+=batch[i];
		          if(_currentMax == batch[i]){
		        	  _currentMaxPos = _currentPos;
		        	  cent=i;
		          }
		        }else{
		          _aboveThreshold = false;
		          _currentMean/=(i-beg+1);
//		          doWrite(batchP, batchM, cent);
		          doWrite(batchP, batchM, (beg+i)/2);
		        }
		      }
		    }
	  }
	  
	  private void doWrite() throws IOException {
		long centerPos = (long)_currentMaxPos - _startPeakPos;
	    bw.write(chr + "\t" + _startPeakPos + "\t" + (_currentPos-1) + "\t" + (chr + "." + _counter++) + "\t" + "0" + "\t" + "." + "\t" + nf.format(_currentMax) + "\t" + "-1" + "\t" + "-1" + "\t" + centerPos + "\n");
	    _currentMax = 0.0f;
	    _currentMaxPos = 0;
	    _startPeakPos = 0l;
	  }
	  
	  private void doWrite(float[] batchP, float[] batchM, int cent) throws IOException {
		  	int win=100, binN=20, inc;
		  	long centerPos = (long)_currentMaxPos - _startPeakPos;
		  	String pmstr="";
		  	double dist2uni, sump, summ;
		  	double[] psig=new double[binN], msig=new double[binN];
		  	inc=2*win/binN;
		  	sump=summ=0.01;
		  	dist2uni=0;
		  	for (int i=0;i<binN;i++){
		  		int beg=cent-win+i*inc;
		  		float meanP,meanM;
		  		meanP=meanM=0;
		  		for (int j=0;j<inc;j++){
			  		int ind=Math.max(0, beg+j);
			  		ind=Math.min(batchP.length-1, ind);
			  		meanP+=batchP[ind];
			  		meanM+=batchM[ind];
		  		}
			  	pmstr+="\t"+Math.log(meanP/meanM);
			  	psig[i]=meanP;
			  	msig[i]=meanM;
			  	sump+=meanP;
			  	summ+=meanM;
		  	}
		  	for (int i=0;i<binN;i++){
		  		double tempP, tempM;
		  		psig[i]/=sump;
		  		msig[i]/=summ;
		  		tempP=psig[i]-1.0/binN;
		  		tempM=msig[i]-1.0/binN;
		  		dist2uni+=(tempP*tempP+tempM*tempM);
		  	}
		  	dist2uni=Math.sqrt(dist2uni)+0.001;
//		    bw.write(chr + "\t" + _startPeakPos + "\t" + (_currentPos-1) + "\t" + (chr + "." + _counter++) + "\t" + "0" + "\t" + "." + "\t" + nf.format(_currentMax) + "\t" + "-1" + "\t" + "-1" + "\t" + centerPos + pmstr+ "\n");
		  	bw.write(chr + "\t" + _startPeakPos + "\t" + (_currentPos-1) + "\t" + (chr + "." + _counter++) + "\t" + nf.format(_currentMean/dist2uni) + "\t" + "." + "\t" + nf.format(_currentMean) + "\t" + "-1" + "\t" + "-1" + "\t" + centerPos +"\n");
		    _currentMax = 0.0f;
		    _currentMaxPos = 0;
		    _startPeakPos = 0l;
	  }
}
