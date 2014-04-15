package io;

import java.io.IOException;

public interface DensityWriter {
  public void writeDensity(float[] batch, int start, int offset) throws IOException;
  
  public void writeDensityPM(float[] batch, int start, int length, float[] batchP, float[] batchM) throws IOException;
  
  public void close() throws IOException;
  
  public void setThreshold(float threshold);
}
