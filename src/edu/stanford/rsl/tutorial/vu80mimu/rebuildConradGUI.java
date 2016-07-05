package edu.stanford.rsl.tutorial.vu80mimu;
import java.io.IOException;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.filtering.ImageFilteringTool;
import edu.stanford.rsl.conrad.filtering.redundancy.TrajectoryParkerWeightingTool;
import edu.stanford.rsl.conrad.io.DennerleinProjectionSource;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.FileInfo;

public class rebuildConradGUI {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String filename = "/proj/i5fpctr/data/Exercise5/8SDR_HORIZONTAL_KNEE_0002.bin";
		DennerleinProjectionSource source = new DennerleinProjectionSource();
		FileInfo info = source.getHeaderInfo(filename);
		source.initStream(filename);
		Grid3D image3 = new Grid3D(info.width,info.height, info.nImages);
		for(int i = 0;i<image3.getSize()[2];i++){
			image3.setSubGrid(i, source.getNextProjection());
		}
		

		
		
		
		
		new ImageJ();
		Configuration conf = Configuration.loadConfiguration("configKnee.xml");	
		
		
		TrajectoryParkerWeightingTool parker = new TrajectoryParkerWeightingTool();
		
		
		ImageFilteringTool filters[] = new ImageFilteringTool[]{parker};
		Grid3D filteredIm = ImageUtil.applyFiltersInParallel(image3, filters);
		filteredIm.show("Parker Only");

	}

}
