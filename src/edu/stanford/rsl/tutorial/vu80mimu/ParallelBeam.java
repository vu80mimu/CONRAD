package edu.stanford.rsl.tutorial.vu80mimu;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Transform;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.FFTUtil;
import edu.stanford.rsl.tutorial.phantoms.SheppLogan;
import ij.ImageJ;

public class ParallelBeam {
	public static Grid2D createSinogram(Grid2D image, int numProjections, double[] spacingDetector, int numDetectorPixels) {
		Grid2D sinogramm = new Grid2D(numProjections,numDetectorPixels); // sino[1] = s und sino[0] = theta
		sinogramm.setSpacing(spacingDetector);
		sinogramm.setOrigin(-(sinogramm.getSize()[0]*spacingDetector[0]/2),-( sinogramm.getSize()[1]*spacingDetector[1]/2));
		// Set sampling rate
		//TODO
		final double samplingStepSize = 0.3;
		// create box around image
		Box box = new Box((image.getSize()[0] * image.getSpacing()[0]), (image.getSize()[1] * image.getSpacing()[1]),
				1);
		Translation transform = new Translation(new double[] { -(image.getSpacing()[0] * image.getSize()[0]) / 2,
				-(image.getSpacing()[1] * image.getSize()[1] )/ 2, -1 });
		box.applyTransform(transform);

		double splitAngle = (180.0/numProjections) *2*Math.PI / 360;
		// go over all projections in 180 degree => column by column
		for(int i = 0; i <numProjections; i++){ // go over the columns
			// calculate angle theta of of the ith projection
			double theta = i * splitAngle;
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);
			// go over single rays
			for(int j = 0 ; j < numDetectorPixels; j++){
				// calculate s in World coordinates
				double s = j* spacingDetector[0] - (numDetectorPixels*spacingDetector[0])/2;
				// define two points on the line through the box
				PointND p1 = new PointND(s * cosTheta, s * sinTheta, .0d); // 
				PointND p2 = new PointND(-sinTheta + (s * cosTheta),
						(s * sinTheta) + cosTheta, .0d);
				
				// set up line equation
				StraightLine line = new StraightLine(p1, p2);
				// compute intersections between bounding box and intersection line.
				ArrayList<PointND> points = box.intersect(line); 
				// only if we have intersections
				if (2 != points.size()){
					if(points.size() == 0) {
						line.getDirection().multiplyBy(-1.d);
						points = box.intersect(line);
					}
					if(points.size() == 0)
						//System.out.println("No points found");
						continue;
				}
				// intersection points
				PointND start = points.get(0); // [mm]
				PointND end = points.get(1);   // [mm]
				//System.out.println("["+ start +"," +end+"]");
				
				// calculate length of intersection line
				double length = Math.sqrt((double) (Math.pow(start.get(0)- end.get(0), 2)+ Math.pow(start.get(1)-end.get(1), 2)));
				//System.out.println(length);
				//intersection points
				int samplingRate = (int)(length/samplingStepSize);
				PointND currentPoint = new PointND(start);
				SimpleVector step = new SimpleVector(line.getDirection().multipliedBy(samplingStepSize));
				// subtract step the first time so that first ray is at the start
				currentPoint.getAbstractVector().subtract(step);
				
				// initialize sum for integral calculation
				float sum = 0.f;
				for(int l = 0; l <= samplingRate; l++){
					currentPoint.getAbstractVector().add(step);
					//double[] currentIndex = image.physicalToIndex(currentPoint.get(0), currentPoint.get(1));
					double currentX = (currentPoint.get(0) / image.getSpacing()[0]) - image.getOrigin()[0];
					double currentY = (currentPoint.get(1) / image.getSpacing()[1]) - image.getOrigin()[1];
					float value = InterpolationOperators.interpolateLinear(image, currentX, currentY);
					//System.out.println(value);
					sum += value;
					//System.out.println(sum);
					
				}
				sinogramm.setAtIndex(i, j, sum); 
				
			}
			
		}
		return sinogramm;
	}
	
	public static Grid2D backProject(Grid2D sino){
		
		if (sino.getSize()[0] < 180){
			System.err.println("This is not a full scan.Please use at least 180°.");
			return null;
		}
		Grid2D image = new Grid2D(sino.getSize()[1],sino.getSize()[1]);
		image.setSpacing(new double []{sino.getSpacing()[1], sino.getSpacing()[1]});
		image.setOrigin(-(image.getSize()[0]*image.getSpacing()[0]/2),-( image.getSize()[1]*image.getSpacing()[1]/2));
		
		for (int t = 0; t < sino.getSize()[0]; t++){
			// calc the actual theta values in deg
			double theta = t*(sino.getSize()[0]/180);
			double cosTheta = Math.cos(theta *(2*Math.PI/360));
			double sinTheta = Math.sin(theta *(2*Math.PI/360));
			for (int i = 0; i < image.getSize()[0]; i++){
				for(int j = 0; j < image.getSize()[1]; j++){
					//pixels to worldcoordinates
					double [] physIndex = (image.indexToPhysical(i, j));
					double s = physIndex[1]*sinTheta + physIndex[0]*cosTheta;
					double [] sinoIndex = sino.physicalToIndex(t,s);
					float val = InterpolationOperators.interpolateLinear(sino, t, sinoIndex[1]);
					float pixelVal = image.getAtIndex(i, j)+ val;
					image.setAtIndex(i, j, pixelVal);
				}
			}
		}
		return image;
	}
	public static Grid2D rampFilter(Grid2D sino){
		Grid2D filtSino = new Grid2D(sino.getSize()[0],sino.getSize()[1]);
		filtSino.setSpacing(new double []{sino.getSpacing()[0], sino.getSpacing()[1]});
		filtSino.setOrigin(-(filtSino.getSize()[0]*filtSino.getSpacing()[0]/2),-( filtSino.getSize()[1]*filtSino.getSpacing()[1]/2));
		// calculate frequency spacing
		Grid1D ramp = new Grid1D(FFTUtil.getNextPowerOfTwo(sino.getSize()[1]));
		int k = ramp.getSize()[0];
		double freqSpacing = 1/(sino.getSpacing()[1]*k);
		for (int i = 0; i< k/2; i++){
			ramp.setAtIndex(i, (float)Math.abs(2.0f*Math.PI*freqSpacing*i));
			ramp.setAtIndex(k-i-1, (float)Math.abs(2.0f*Math.PI*freqSpacing*i));
		}
		ramp.show();
		Grid1DComplex complexRamp = new Grid1DComplex(ramp);
		// FFT of one line of the sino
		
		for(int i = 0; i < sino.getSize()[0]; i++){
			Grid1DComplex complexColumn = new Grid1DComplex(sino.getSize()[1]);
			for(int j = 0; j < sino.getSize()[1]; j++){
				complexColumn.setAtIndex(j, sino.getAtIndex(i, j));
			}
			complexColumn.transformForward();

			// multiply with ramp filter
			for(int j = 0 ; j < complexColumn.getSize()[0]; j++){
				complexColumn.multiplyAtIndex(j, complexRamp.getRealAtIndex(j), complexRamp.getImagAtIndex(j));
			}
			
			complexColumn.transformInverse();
			// write filtered sinogramm into output image
			for(int j = 0 ; j < sino.getSize()[1]; j++){
				filtSino.setAtIndex(i, j, complexColumn.getRealAtIndex(j));
			}
		}
		
		
	return filtSino;	
	}
	
	public static Grid2D ramLakFilter(Grid2D sino){
		Grid2D filtSino = new Grid2D(sino.getSize()[0],sino.getSize()[1]);
		filtSino.setSpacing(new double []{sino.getSpacing()[0], sino.getSpacing()[1]});
		filtSino.setOrigin(-(filtSino.getSize()[0]*filtSino.getSpacing()[0]/2),-( filtSino.getSize()[1]*filtSino.getSpacing()[1]/2));
		// calculate frequency spacing
		Grid1D ramLak = new Grid1D(sino.getSize()[1]);
		Grid1DComplex complexRamLak = new Grid1DComplex(ramLak);
		int k = complexRamLak.getSize()[0];
		//double freqSpacing = 1/(sino.getSpacing()[0]*k);
		for (int i = 0; i<= k/2; i++){
			if(i == 0){
				complexRamLak.setAtIndex(i, 0.25f);
			}else if(i%2 == 0){
				// even
				complexRamLak.setAtIndex(i, 0.0f);
				complexRamLak.setAtIndex(k-i, 0.0f);
			}else{
				complexRamLak.setAtIndex(i, (float)(-1/(Math.pow(i, 2)*Math.pow(Math.PI, 2))));
				complexRamLak.setAtIndex(k-i, (float)(-1/(Math.pow(i, 2)*Math.pow(Math.PI, 2))));
			}

		}
		complexRamLak.show();
		//Grid1DComplex complexRamLak = new Grid1DComplex(ramLak);
		complexRamLak.transformForward();
		complexRamLak.show();
		
		// FFT of one line of the sino
		
		for(int i = 0; i < sino.getSize()[0]; i++){
			Grid1DComplex complexColumn = new Grid1DComplex(sino.getSize()[1]);
			for(int j = 0; j < sino.getSize()[1]; j++){
				complexColumn.setAtIndex(j, sino.getAtIndex(i, j));
			}
			complexColumn.transformForward();

			// multiply with ramp filter
			for(int j = 0 ; j < complexColumn.getSize()[0]; j++){
				complexColumn.multiplyAtIndex(j, complexRamLak.getRealAtIndex(j), complexRamLak.getImagAtIndex(j));
			}
			
			complexColumn.transformInverse();
			// write filtered sinogramm into output image
			for(int j = 0 ; j < sino.getSize()[1]; j++){
				filtSino.setAtIndex(i, j, complexColumn.getRealAtIndex(j));
			}
		}
		
		
	return filtSino;	
	}
	public static Grid2D filteredBackprojection(Grid2D sino, String filterType){
		Grid2D imOut = new Grid2D(sino.getSize()[0],sino.getSize()[1]);
		imOut.setSpacing(new double []{sino.getSpacing()[0], sino.getSpacing()[1]});
		imOut.setOrigin(-(imOut.getSize()[0]*imOut.getSpacing()[0]/2),-( imOut.getSize()[1]*sino.getSpacing()[1]/2));
		switch(filterType){
		case "ramp":
			imOut = rampFilter(sino);
			break;
		case "ramLak":
			imOut = ramLakFilter(sino);
			break;
		default:
			System.err.println("The filter type " + filterType +" is not implemented.");
		}
		
		imOut = backProject(imOut);
		
		return imOut;
	}
	public static void main(String args[]) {
		new ImageJ();
		// 1. Create the Shepp Logan Phantom
		//SheppLogan sheppLoganPhantom = new SheppLogan(256);
		//sheppLoganPhantom.show();
		//Create Phantom
		CustPhantom phantom = new CustPhantom(200,300 , new double[] { 1.0, 1.0 });
		//ParallelBeam p = new ParallelBeam();
		phantom.show();
		//Create Sinogramm of phantom
		Grid2D sinogramm = createSinogram(phantom, 180, new double[] { 1.0, 1.0 }, 366);
		sinogramm.show("Sinogramm");
		//perform backprojection without filtering
		//Grid2D backProj = backProject(sinogramm);
		//backProj.show("Back Projection");
		
		Grid2D rampFil = rampFilter(sinogramm);
		rampFil.show("Ramp Filt");
		
		Grid2D ramlakFil = ramLakFilter(sinogramm);
		ramlakFil.show("Ram Lak");
		
		//perform filteredBackprojection
		Grid2D filtBackProj = filteredBackprojection(sinogramm, "ramp");
		filtBackProj.show("Filtered Backprojection");
		
		Grid2D filtBackProjR = filteredBackprojection(sinogramm, "ramLak");
		filtBackProjR.show("Filtered Backprojection with RamLak");
	}

}
