// Texture sampling
__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
__kernel void backproject(int sizeSino, int sizeImage, __read_only image2d_t sinogram, __global float *image){
	
	
	int indSinox = get_global_id(0);
	int indSinoy = get_gloaal_id(1);
	
	int indImage = get_global_id(0);
	
	
	if(indSino >= sizeSino) return;
	if(indImage >= sizeImage) return;

	
	//phan[ind]= phan[ind] + phan2[ind];

}