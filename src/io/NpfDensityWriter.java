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
  private final int offset;
  private final String chr;
  
  private float _threshold = 0.0f;
  private long _currentPos;
  private long _counter = 1;
  
  private long _startPeakPos = 0;
  private boolean _aboveThreshold = false;
  private float _currentMax = 0.0f;
  private long _currentMaxPos = 0;
  private float _currentMean = 0.0f;
  private float _currentP = 0.0f;
  private float _currentM = 0.0f;
  
  private float[] _queP, _queM;
  
  private float _currentMaxL=0, _currentMaxR=0;
  private long _currentPosL=0, _currentPosR=0;
  private double _currentRatioL = 0;
  private double _currentRatioR = 0;
  
  private NumberFormat nf;
  
  public NpfDensityWriter(File f, String chr, long chromStart, int offset) throws IOException {
	    bw = new BufferedWriter(new FileWriter(f));
	    this.chromPos = chromStart;
	    this.chr = chr;
	    this.offset = offset;
	    nf = NumberFormat.getNumberInstance();
	    nf.setGroupingUsed(false);
	    nf.setMaximumFractionDigits(8);
	    nf.setMinimumFractionDigits(8);
	    
	    _queP=new float[offset];
	    _queM=new float[offset];
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
	  
	  public int getMaxLeft(float[] batch,float[] batchP, float[] batchM, int cent, int offset){
		  int head, tail, maxIdx=cent-offset/2;
		  double sumP=0.001, sumM=0.001, maxR=1.0;
		  tail=cent;
		  for (head=cent; head>=0; head--){
			  if (batch[head]<_threshold/2)
				  break;
			  sumP+=batchP[head];
			  sumM+=batchM[head];
			  if (tail-head>=offset){
				  if (sumP/sumM>maxR){
					  maxR=sumP/sumM;
					  maxIdx=(head+tail)/2;
				  }
				  sumP-=batchP[tail];
				  sumM-=batchM[tail];
				  tail--;
			  }
		  }
		  _currentRatioL=maxR;
		  return maxIdx;
	  }
	  
	  public int getMaxRight(float[] batch,float[] batchP, float[] batchM, int cent, int offset){
		  int head, tail, maxIdx=cent+offset/2;
		  double sumP=0.001, sumM=0.001, maxR=1.0;
		  tail=cent;
		  for (head=cent; head<batch.length; head++){
			  if (batch[head]<_threshold/2)
				  break;
			  sumP+=batchP[head];
			  sumM+=batchM[head];
			  if (head-tail>=offset){
				  if (sumM/sumP>maxR){
					  maxR=sumM/sumP;
					  maxIdx=(head+tail)/2;
				  }
				  sumP-=batchP[tail];
				  sumM-=batchM[tail];
				  tail++;
			  }
		  }
		  _currentRatioR=maxR;
		  return maxIdx;
	  }
	  
	  public void writeDensityPM(float[] batch, int start, int length, float[] batchP, float[] batchM)
		      throws IOException {
		    int end = start + length;
		    int cent=0,beg=0, tail;
		    float tempRatio;
		    long left, right;
		    long batchStart=_currentPos;
		    for(int i = start; i < end; ++i, ++_currentPos){
		      if(!_aboveThreshold){
		        if(batch[i] > _threshold){
		          _aboveThreshold = true;
		          _startPeakPos = _currentPos;
		          _currentMax = batch[i];
		          _currentMaxPos = _currentPos;
		          _currentMean = batch[i];
		          tail=(int)(_currentPos-_startPeakPos)%offset;
		          _queP[tail]=batchP[i];
		          _queM[tail]=batchM[i];
		          _currentP=0.001f+batchP[i];
		          _currentM=0.001f+batchM[i];
		          _currentMaxL=1.0f;
		          _currentMaxR=1.0f;
		          _currentPosL=_currentPosR=_currentPos;
		          cent=i;
		          beg=i;
		        }
		      }else{ // aboveThreshold
		        if(batch[i] > _threshold){
		          _currentMax = Math.max(_currentMax, batch[i]);
		          _currentMean+=batch[i];
		          tail=(int)(_currentPos-_startPeakPos)%offset;
		          if (_currentPos-_startPeakPos>=offset){
		        	  tempRatio=_currentP/_currentM;
		        	  if (tempRatio>_currentMaxL){
		        		  _currentMaxL=tempRatio;
		        		  _currentPosL=_currentPos-offset/2;
		        	  }
		        	  if (tempRatio<_currentMaxR){
		        		  _currentMaxR=tempRatio;
		        		  _currentPosR=_currentPos-offset/2;
		        	  }
		        	  _currentP-=_queP[tail];
		        	  _currentM-=_queM[tail];
		          }
		          _queP[tail]=batchP[i];
		          _queM[tail]=batchM[i];
		          _currentP+=batchP[i];
		          _currentM+=batchM[i];
		          if(_currentMax == batch[i]){
		        	  _currentMaxPos = _currentPos;
		        	  cent=i;
		          }
		        }else{
		          _aboveThreshold = false;
		          _currentMean/=(_currentPos-_startPeakPos+1);
//		          _currentP/=(i-beg+1);
//		          _currentM/=(i-beg+1);
		          // issue: may have overlap
//		          left=getMaxLeft(batch, batchP, batchM, (beg+i)/2, offset);
//		          right=getMaxRight(batch, batchP, batchM, (beg+i)/2, offset);
//		          _currentP=0;
//		          _currentM=0;
//		          for (int j=(int)left;j<=(int)right;j++){
//		        	  if (j<0||j>=batchP.length)
//		        		  continue;
//		        	  _currentP+=batchP[j];
//		        	  _currentM+=batchM[j];
//		          }
//		          left+=batchStart;
//		          right+=(batchStart+1);
		          if (_currentPosL<_currentPosR){
		        	  left=_currentPosL;
		        	  right=_currentPosR;
		          }
		          else{
		        	  _currentMaxL=_currentMaxR=1.0f;
		        	  left=_startPeakPos;
		        	  right=_currentPos;
		          }
		          doWrite(left, right);
//		          doWrite(Math.min(_startPeakPos, left), right);
//		          doWrite(batchP, batchM, (beg+i)/2);
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
	  
	  private void doWrite(long left, long right) throws IOException {
			long centerPos = (long)_currentMaxPos - left;
			// issue: _currentMaxPos may not be in the region
		    bw.write(chr + "\t" + left + "\t" + (right-1) + "\t" + (chr + "." + _counter++) + "\t" + nf.format(_currentMaxL) + "\t" + nf.format(1/_currentMaxR) + "\t" + nf.format(_currentMean) + "\t" + "-1" + "\t" + "-1" + "\t" + centerPos + "\n");
//			bw.write(chr + "\t" + left + "\t" + (right-1) + "\t" + (chr + "." + _counter++) + "\t" + nf.format(_currentP) + "\t" + nf.format(_currentM) + "\t" + nf.format(_currentMax) + "\t" + "-1" + "\t" + "-1" + "\t" + centerPos + "\n");
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
