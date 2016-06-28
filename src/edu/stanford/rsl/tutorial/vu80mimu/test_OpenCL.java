package edu.stanford.rsl.tutorial.vu80mimu;
import ij.ImageJ;

import java.io.IOException;
import java.nio.FloatBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLMemory.Mem;
import com.jogamp.opencl.CLProgram;
import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;


public class test_OpenCL {
	
	public static Grid2D addGPU(int width, int height){
		
		CustPhantom phantom = new CustPhantom(width,height,new double[]{1.0,1.0});
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice clDevice = context.getMaxFlopsDevice();
		
		float[] buffer1 = phantom.getBuffer().clone(); // nicht auf phanom arbeiten sondern buffer neu erstellen
		CLBuffer<FloatBuffer> clPhantomFloat = context.createFloatBuffer(buffer1.length,Mem.READ_WRITE);
		clPhantomFloat.getBuffer().put(buffer1);
		clPhantomFloat.getBuffer().rewind();
		
		//CustPhantom phantom2 = new CustPhantom(width,height,new double[]{1.0,1.0});
		
		float[] buffer2 = phantom.getBuffer().clone();
		
		CLBuffer<FloatBuffer> clPhantomFloat2 = context.createFloatBuffer(buffer2.length, Mem.READ_WRITE);
		clPhantomFloat2.getBuffer().put(buffer2);
		clPhantomFloat2.getBuffer().rewind();
		
		// load and complite cl-code that contains kernel methods
		CLProgram program = null;
		try{
			program = context.createProgram(test_OpenCL.class.getResourceAsStream("addPhantomsCL.cl")).build();
		}catch(IOException ex){
			ex.printStackTrace();
		}
		// kernel function with input parameter
		CLKernel kernelFunction = program.createCLKernel("add");
		kernelFunction.rewind(); //pointer geht an den Anfang
		
		
		int localWorkSize = 128;
		int globalWorkSize = OpenCLUtil.roundUp(localWorkSize, buffer1.length);
		
		CLCommandQueue queue = clDevice.createCommandQueue();
		// writing into phantom1
		queue.putWriteBuffer(clPhantomFloat, true).finish();
		
		kernelFunction.putArg(width*height)
		.putArg(clPhantomFloat).putArg(clPhantomFloat2); // in addPhantomCl OpenCLGrid2D ??
		
		queue.put1DRangeKernel(kernelFunction, 0, globalWorkSize, localWorkSize)
		.finish(); // fuehre kernel function aus
		
		queue.putReadBuffer(clPhantomFloat, true);
		
		for(int j = 0; j< phantom.getSize()[1]; j++){
		for(int i = 0; i< phantom.getSize()[0]; i++){
			phantom.setAtIndex(i, j, clPhantomFloat.getBuffer().get()); // buffer reads column wise
		}
		}
		
		return phantom;
		
	}
	
public static Grid2D backProjectGPU(Grid2D sino, int localWorkSize){
		
		if (sino.getSize()[0] < 180){
			System.err.println("This is not a full scan.Please use at least 180°.");
			return null;
		}
		
		
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice clDevice = context.getMaxFlopsDevice();
		
		
		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY,ChannelType.FLOAT);
		OpenCLGrid2D sinoCl = new OpenCLGrid2D(sino);
		CLImage2d<FloatBuffer> sinoTexture = null;
		sinoTexture = context.createImage2d(sinoCl.getDelegate().getCLBuffer().getBuffer(), sinoCl.getSize()[0], sinoCl.getSize()[1], format,Mem.READ_ONLY );
		
				
		
		Grid2D image = new Grid2D(sino.getSize()[1],sino.getSize()[1]);
		image.setSpacing(new double []{sino.getSpacing()[1], sino.getSpacing()[1]});
		image.setOrigin(-(image.getSize()[0]*image.getSpacing()[0]/2),-( image.getSize()[1]*image.getSpacing()[1]/2));
		
		
		
		float[] bufferImage = image.getBuffer().clone(); // nicht auf phanom arbeiten sondern buffer neu erstellen
		CLBuffer<FloatBuffer> clImageFloat = context.createFloatBuffer(bufferImage.length,Mem.READ_WRITE);
		clImageFloat.getBuffer().put(bufferImage);
		clImageFloat.getBuffer().rewind();
		
		// load and complite cl-code that contains kernel methods
		CLProgram program = null;
		try{
			program = context.createProgram(test_OpenCL.class.getResourceAsStream("backProject.cl")).build();
		}catch(IOException ex){
			ex.printStackTrace();
		}
		// kernel function with input parameter
		CLKernel kernelFunction = program.createCLKernel("backproject");
		kernelFunction.rewind(); //pointer geht an den Anfang
		
		int globalWorkSize = OpenCLUtil.roundUp(localWorkSize, bufferImage.length); // sample in output space		
		
		
		CLCommandQueue queue = clDevice.createCommandQueue();
		// memory transfer
		queue.putWriteImage(sinoTexture, true).putWriteBuffer(clImageFloat,true).finish();
		
		kernelFunction.putArg(image.getSize()[0]).putArg(image.getSize()[1])
		.putArg(sino.getSize()[0]).putArg(sino.getSize()[1]).putArg(sinoTexture).putArg(clImageFloat)
		.putArg((int)image.getOrigin()[0]).putArg((int)image.getOrigin()[1]).putArg((float) image.getSpacing()[0])
		.putArg((float)image.getSpacing()[1]).putArg((int) sino.getOrigin()[1]).putArg((float) sino.getSpacing()[1]); // in addPhantomCl OpenCLGrid2D ??
		
		queue.put1DRangeKernel(kernelFunction, 0, globalWorkSize, localWorkSize)
		.finish(); // fuehre kernel function aus
		
		queue.putReadBuffer(clImageFloat, true);
		
		for(int j = 0; j< image.getSize()[1]; j++){
		for(int i = 0; i< image.getSize()[0]; i++){
			image.setAtIndex(i, j, clImageFloat.getBuffer().get()); // buffer reads column wise
		}
		}
		
		return image;
	}
	
	

	public static void main(String[] args) {
		/*
		// TODO Auto-generated method stub
		CustPhantom phantom = new CustPhantom(256,256,new double[]{1.0,1.0});
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice clDevice = context.getMaxFlopsDevice();
		OpenCLGrid2D phantom2 = new OpenCLGrid2D(phantom, context, clDevice);
		//OpenCLGrid2D phantom2 = new OpenCLGrid2D(phantom);
		int num = 1;
		
		long startCPU = System.currentTimeMillis();
		for(int i = 0; i<num; i++){
			NumericPointwiseOperators.addBy(phantom, phantom);
			
		}
		long endCPU = System.currentTimeMillis();
		
		long startGPU = System.currentTimeMillis();
		for(int i = 0; i<num; i++){
			NumericPointwiseOperators.addBy(phantom2, phantom2);
			
		}
		long endGPU = System.currentTimeMillis();
		
		System.out.println("CPU time:" + (endCPU -startCPU));
		System.out.println("GPU time:" + (endGPU -startGPU));
		
		new ImageJ();
		Grid2D phan = addGPU(256,256);
		phan.show();
		
		phantom.show();
		*/
		
		/* Aufgabe 3*/
		
		ParallelBeam pb = new ParallelBeam();
		CustPhantom backprojphantom = new CustPhantom(200,300 , new double[] { 1.0, 1.0 });
		//ParallelBeam p = new ParallelBeam();
		backprojphantom.show();
		//Create Sinogramm of phantom
		Grid2D sinogramm = pb.createSinogram(backprojphantom, 180, new double[] { 1.0, 1.0 }, 366);
		sinogramm.show("Sinogramm");
		int localWorksize = 128;
		Grid2D backproj = backProjectGPU(sinogramm,localWorksize);
		backproj.show();
	}

}
