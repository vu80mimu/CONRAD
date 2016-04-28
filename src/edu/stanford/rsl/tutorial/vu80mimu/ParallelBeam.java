package edu.stanford.rsl.tutorial.vu80mimu;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.tutorial.phantoms.SheppLogan;
import ij.ImageJ;

public class ParallelBeam {
	public Grid2D createSinogram(Grid2D image, int numProjections, double[] spacingDetector, int numDetectorPixels) {
		Grid2D sinogramm = new Grid2D(numDetectorPixels,numProjections);
		// Set sampling rate
		//TODO
		final double samplingStepSize = 3.0;
		// create box around image
		Box box = new Box((image.getSize()[0] * image.getSpacing()[0]), (image.getSize()[1] * image.getSpacing()[1]),
				1);
		Translation transform = new Translation(new double[] { -(image.getSpacing()[0] * image.getSize()[0]) / 2,
				-(image.getSpacing()[1] * image.getSize()[1] )/ 2, -1 });
		box.applyTransform(transform);
		System.out.println("HALLo");
		double splitAngle = 180.0/numProjections;
		// go over all projections in 180 degree
		for(int i = 0; i <numProjections; i++){
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
						System.out.println("No points found");
						continue;
				}
				// intersection points
				PointND start = points.get(0); // [mm]
				PointND end = points.get(1);   // [mm]
				System.out.println("["+ start +"," +end+"]");
				
				// calculate length of intersection line
				double length = Math.sqrt((double) (Math.pow(start.get(0)- end.get(0), 2)+ Math.pow(start.get(1)-end.get(1), 2)));
				System.out.println(length);
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
					double currentX = (currentPoint.get(0) / image.getSpacing()[0]) + 128;
					double currentY = (currentPoint.get(1) / image.getSpacing()[1])+ 128;
					float value = InterpolationOperators.interpolateLinear(image, currentX, currentY);
					//System.out.println(value);
					sum += value;
					//System.out.println(sum);
					
				}
				sinogramm.setAtIndex(j, i, sum);
				
			}
			
		}

		return sinogramm;
	}

	public static void main(String args[]) {
		new ImageJ();
		// 1. Create the Shepp Logan Phantom
		SheppLogan sheppLoganPhantom = new SheppLogan(256);
		sheppLoganPhantom.show();
		CustPhantom phantom = new CustPhantom(256, 256, new double[] { 1.0, 1.0 });
		ParallelBeam p = new ParallelBeam();
		phantom.show();
		Grid2D sinogramm = p.createSinogram(sheppLoganPhantom, 180, new double[] { 1.0, 1.0 }, 256);
		sinogramm.show("Sinogramm");
	}

}
