package edu.stanford.rsl.tutorial.vu80mimu;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import ij.ImageJ;

public class FanBeam {
	public static Grid2D createFanogram(Grid2D image, double[] spacingDetector, int numDetectorPixels, double incBeta, int numProjections, double dSI, double dSD) {
		// ensure the phantom is not hit
		double maxLength = Math.sqrt((Math.pow(image.getSize()[0], 2))+ (Math.pow(image.getSize()[1], 2)));
		if((maxLength/2 >= dSI) || (maxLength/2 >= (dSD - dSI))){
			System.err.println("The space between source and detector is too narrow.");
			return null;
		}
		
		
		Grid2D fanogramm = new Grid2D(numProjections,numDetectorPixels); // sino[1] = s und sino[0] = Beta
		fanogramm.setSpacing(spacingDetector);
		fanogramm.setOrigin(-(fanogramm.getSize()[0]*spacingDetector[0]/2),-( fanogramm.getSize()[1]*spacingDetector[1]/2));
		// Set sampling rate
		//TODO
		final double samplingStepSize = 0.5; // how often are values evaluated along one line Integral
		// create box around image
		Box box = new Box((image.getSize()[0] * image.getSpacing()[0]), (image.getSize()[1] * image.getSpacing()[1]),
				1);
		Translation transform = new Translation(new double[] { -(image.getSpacing()[0] * image.getSize()[0]) / 2,
				-(image.getSpacing()[1] * image.getSize()[1] )/ 2, -1 });
		box.applyTransform(transform);

		double totalRotationAngle = (incBeta * numProjections);// *2*Math.PI / 360;
		// go over all projections in 180 degree => column by column
		for(int i = 0; i <numProjections; i++){ // go over the columns
			// calculate angle Beta of of the ith projection
			double beta = ((i * incBeta))*2*Math.PI / 360;
			double cosBeta = Math.cos(beta);
			double sinBeta = Math.sin(beta);
			// go over single rays
			for(int j = 0 ; j < numDetectorPixels; j++){
				// calculate s in World coordinates
				double t = j* spacingDetector[0] - (numDetectorPixels*spacingDetector[0])/2;
				// define two points on the line through the 3D box
				PointND p1 = new PointND(-sinBeta* dSI, cosBeta*dSI, 0.0d); // 
				PointND midDetector = new PointND(sinBeta*(dSD - dSI), -cosBeta*(dSD -dSI));
				
				PointND p2 = new PointND(midDetector.getCoordinates()[0]+cosBeta*t, midDetector.getCoordinates()[1]+sinBeta*t, 0.0d);
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
				//getDirection returns a non normalized vector
				SimpleVector step = new SimpleVector(line.getDirection().normalizedL2().multipliedBy(samplingStepSize));
				// subtract step the first time so that first ray is at the start
				currentPoint.getAbstractVector().subtract(step);
				
				// initialize sum for integral calculation
				float sum = 0.f;
				for(int l = 0; l <= samplingRate; l++){
					currentPoint.getAbstractVector().add(step);
					//double[] currentIndex = image.physicalToIndex(currentPoint.get(0), currentPoint.get(1));
					double currentX = (currentPoint.get(0) - image.getOrigin()[0]) / image.getSpacing()[0];
					double currentY = (currentPoint.get(1) - image.getOrigin()[1]) / image.getSpacing()[1];
					float value = InterpolationOperators.interpolateLinear(image, currentX, currentY);
					//System.out.println(value);
					sum += value;
					//System.out.println(sum);
					
				}
				fanogramm.setAtIndex(i,j, sum); 
				
			}
			
		}
		return fanogramm;
	}
	
	public static Grid2D rebinning(Grid2D fano, double incBeta,double dSI, double dSD){
		Grid2D sinogramm = new Grid2D(fano.getSize()[0],fano.getSize()[1]); // sino[1] = s und sino[0] = Beta
		sinogramm.setSpacing(fano.getSpacing()[0],fano.getSpacing()[1]);
		sinogramm.setOrigin(-(sinogramm.getSize()[0]*sinogramm.getSpacing()[0]/2),-( sinogramm.getSize()[1]*sinogramm.getSpacing()[1]/2));
		/* theta = gamma + beta
		 * s = dSI sin(gamma)
		 * tan(gamma) = t/ dSI + dSD
		*/
		for(int i = 0; i < sinogramm.getSize()[0]; i++){ //theta
			for(int j = 0; j< sinogramm.getSize()[1]; j++){ // s
				double[] physIndex = sinogramm.indexToPhysical(i, j); // physIndex[1] = s; physIndex[0] = theat
				// gamma = sin^-1(s/DSI) 
				double gamma = Math.asin((double)(physIndex[1]/dSI)); //in rad
				// t = tan(gamma)* (dSI +dSD)
				double t = Math.tan(gamma) * (dSI+dSD);
				// beta = theta - gamma
				double theta = physIndex[0]*180/sinogramm.getSize()[0]; // grad
				gamma = gamma *360/(2 * Math.PI); // grad
				double beta = theta - gamma; // in grad
				double deltaBeta = beta/incBeta;
				double [] index = fano.physicalToIndex(deltaBeta, t);
				float val = InterpolationOperators.interpolateLinear(fano,index[0],index[1]);
				
				sinogramm.setAtIndex(i, j, val);
			}
		}
		
		return sinogramm;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		new ImageJ();
		CustPhantom phantom = new CustPhantom(200,300 , new double[] { 1.0, 1.0 });
		Grid2D fano = createFanogram(phantom, new double[] {1.0,1.0}, 600, 1.0, 160, 500,1000 );
		phantom.show();
		fano.show("Fanogramm");
		
		
		Grid2D sino = rebinning(fano, 1.0, 500, 1000);
		sino.show("Sinogramm");
		
		ParallelBeam p = new ParallelBeam();
		Grid2D fbp = p.filteredBackprojection(sino, "ramLak");
		fbp.show("Filtered Backproject RamLak");

	}
	
	
	

}
