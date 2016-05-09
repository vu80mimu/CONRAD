package edu.stanford.rsl.tutorial.vu80mimu;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.NumericGrid;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;

public class CustPhantom extends Grid2D{

	public CustPhantom(int width, int height, double[] spacing) {
		super(width, height);
		
		
		
		
		// TODO Auto-generated constructor stub
		if (width >= height){
			this.createEllipse((int)((width*0.4)),(int)(height*0.3), (int)(width/2), (int)(height/2),(float)0.8);
			this.createRect((int)(width*0.1), (int)(height*0.4), (int)(width/2-(width*0.1)/2), (int)(height*0.3),(float) 0.2 );
			this.createCircle(10, width/4, height/2, (float)1);
			this.createCircle(10, width-width/4, height/2, (float)1);
			this.createCircle((int)(3*height/16),(int)(3*width/4), (int) (3*height/4),  (float)0.5);
		}else{
			this.createEllipse((int)((width*0.3)),(int)(height*0.4), (int)(width/2), (int)(height/2),(float)0.8);
			this.createRect((int)(width*0.1), (int)(height*0.4), (int)(width/2-(width*0.1)/2), (int)((height/2)-((height*0.4)/2)),(float) 0.2 );
			this.createCircle(10, width/4, height/2, (float)1);
			this.createCircle(10, width-width/4, height/2, (float)1);
			this.createCircle((int)(3*width/16),(int)(3*width/4), (int) (height/4),  (float)0.5);
		}
		
		this.setSpacing(spacing);
		// changed to minus
		this.setOrigin(-(this.getSize()[0]*spacing[0]/2), -(this.getSize()[1]*spacing[1]/2));

	}
	public void createEllipse(int a,int b,int centerX, int centerY, float value){

		
		int originX = centerX - a;
		int originY = centerY - b;
		for(int i = 0; i < 2*a; i++){
			for(int j = 0; j < 2*b; j++){

				if(Math.pow((double)(i-a),2)/Math.pow((double) a, 2)+Math.pow((double)(j-b),2)/Math.pow((double)b, 2) < 1){
					if(this.getAtIndex(originX+i, originY+j) != (float)0.0){
						float adjustedIntensity = (this.getAtIndex(originX+i, originY+j) + value)/2;
						this.setAtIndex(originX+i, originY+j, adjustedIntensity);
					}else{
						this.setAtIndex(originX+i, originY+j, value);
					}
					
				}
			}
		}
	}
	public void createRect(int w, int h, int originX, int originY, float value){
		for(int j = 0; j < w; j++){
			for(int i = 0; i < h; i++){
				if(this.getAtIndex(originX+j, originY+i) != (float) 0.0){
					float adjustedIntensity = (this.getAtIndex(originX+j, originY+i) + value)/2;
					this.setAtIndex(j+originX,i+originY, adjustedIntensity);
				}else{
					this.setAtIndex(j+originX,i+originY, value);
				}

			}
		}
	}
	public void createCircle(int r, int centerX, int centerY, float value){
		assert centerX > r : "Ellipse Size exeeds image dimensions. Please adjust centerX and/or Radius";
		assert centerY > r : "Ellipse Size exeeds image dimensions. Please adjust centerX and/or Radius";
		int originX = centerX - r;
		int originY = centerY - r;
		for(int i = 0; i < 2*r; i++){
			for(int j = 0; j < 2*r; j++){
				if((Math.pow((double)(i-r),2)+Math.pow((double)(j-r),2)) < Math.pow((double)r, 2)){
					if(this.getAtIndex(originX+i, originY+j) != (float)0.0){
						float adjustedIntensity = (this.getAtIndex(originX+i, originY+j) + value)/2;
						this.setAtIndex(originX+i, originY+j, adjustedIntensity);
					}else{
						this.setAtIndex(originX+i, originY+j, value);
					}

				}
			}
		}
	}
	
	public static void main(String args[]){
		new ImageJ();
		
		CustPhantom phantom = new CustPhantom(256,256,new double[]{1.0,1.0});

		Grid2D image = new Grid2D(256,256);
		for(int j = 0; j < 50; j++){
			for(int i = 0; i < 25; i++){
				image.setAtIndex(j+image.getWidth()/2,i+ image.getHeight()/2, (float)0.1);
			}
		}
		
		
		//Grid2D add = new Grid2D(image);
		//NumericPointwiseOperators.addedBy(add,phantom);
		
		NumericGrid add = NumericPointwiseOperators.addedBy(image, phantom);

		image.show("Bar");
		phantom.show("Phantom");
		add.show("Addition");
	
	
		CustPhantom phantom2 = new CustPhantom(200,256, new double[]{1.0,1.0});
		CustPhantom phantom3 = new CustPhantom(256,200, new double[]{1.0,1.0});
	
		phantom2.show("Phantom2");
		phantom3.show("Phantom3");
		
		NumericGrid add2 = NumericPointwiseOperators.addedBy(phantom3,image );
		add2.show("Addition2");
	
	}

}


//Grid2D phantom = new Grid2D(width, height);
