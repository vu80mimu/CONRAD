// Texture sampling
__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
__kernel void backproject(int sizeImageX,int sizeImageY, int sizeSinoX, int sizeSinoY, __read_only image2d_t sinogram, __global float *image, int originX, int originY, float spacingX, float spacingY, int sinoOriginY, float sinoSpacingY ){
	
	
// 	int indSinox = get_global_id(0);
// 	int indSinoy = get_global_id(1);
	
	int indImagex = get_global_id(0);
	
/*	
	if(indSinox >= sizeSinoX) return;
	if(indSinoy >= sizeSinoY) return;*/
	
	if(indImagex >= sizeImageX) return;
	//if(indImagey >= sizeImageY) return;

	

	
	
	
	
		float valSum = 0;
		for (int t = 0; t < sizeSinoX; t++){
			// calc the actual theta values in deg
			float theta = t*(sizeSinoX/180);
			float cosTheta = cos(theta *(2*M_PI_F/360));
			float sinTheta = sin(theta *(2*M_PI_F/360));
			//pixels to worldcoordinates
			
			
// 			//double [] physIndex = (image.indexToPhysical(i, j));
			int i = indImagex / sizeImageY;
			int j = indImagex % sizeImageY;
			
			float physIndX = i * spacingX + originX;
 			float physIndY = j * spacingY + originY;
			
			float s = physIndY * sinTheta + physIndX*cosTheta;
			
			int sinoIndexY = (s - sinoOriginY)/sinoSpacingY;
			float val = read_imagef(sinogram, sampler, (float2) (t+0.5f, sinoIndexY+0.5f)).x;
// 			float val = InterpolationOperators.interpolateLinear(sino, t, sinoIndex[1]);
// 			float pixelVal = image.getAtIndex(i, j)+ val;
 			valSum += val;
			
			//image[0] = image[0] +1; 
		}
		//image.setAtIndex(i, j, pixelVal);
		image[indImagex] = valSum;

}