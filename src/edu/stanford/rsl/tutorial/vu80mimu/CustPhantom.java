package edu.stanford.rsl.tutorial.vu80mimu;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;

public class CustPhantom extends Grid2D{

	public CustPhantom(int width, int height) {
		super(width, height);
		
		
		
		
		// TODO Auto-generated constructor stub
		
		
		this.createRect((int)(width/4), (int)(height/8), (int)(width/2), (int)(height/16) );
		this.createCircle(50, 50, 50);
		this.createEllipse(50,20 , 100, 100);
		//this.createRect(50, 50, 50, 50);
		this.setAtIndex(50, 50, (float)0.5);
		
	}
	public void createEllipse(int a,int b,int centerX, int centerY){
		int originX = centerX - a;
		int originY = centerY - b;
		for(int i = 0; i < 2*a; i++){
			for(int j = 0; j < 2*b; j++){
				if((Math.pow((double)(i-a),2)/Math.pow(a, 2)+Math.pow((double)(j-b),2))/Math.pow(b, 2) <= 1){
					this.setAtIndex(originX+i, originY+j, (float)0.7);
				}
			}
		}
	}
	public void createRect(int w, int h, int originX, int originY){
		for(int j = originX; j < originX + w; j++){
			for(int i = originY; i < originY + h; i++){
				this.setAtIndex(j,i,(float) 0.3);
			}
		}
	}
	public void createCircle(int r, int centerX, int centerY){
		int originX = centerX - r;
		int originY = centerY - r;
		for(int i = 0; i < 2*r; i++){
			for(int j = 0; j < 2*r; j++){
				if((Math.pow((double)(i-r),2)+Math.pow((double)(j-r),2)) <= Math.pow(r, 2)){
					this.setAtIndex(originX+i, originY+j, (float)0.7);
				}
			}
		}
	}
	public static void main(String args[]){
		new ImageJ();
		CustPhantom phantom = new CustPhantom(256,256);
		phantom.show();
	}
	

}


//Grid2D phantom = new Grid2D(width, height);
