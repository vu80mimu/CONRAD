__kernel void add(int size, __global float *phan, __global float *phan2){
	int ind = get_global_id(0);
	
	
	if(ind >= size) return;

	
	phan[ind]= phan[ind] + phan2[ind];

}