import java.util.ArrayList;
import java.util.Iterator;
import java.util.Vector;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Arrow;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class MotionVectorField_ implements PlugInFilter {

	private ImagePlus imp;
	private int width, height, depth, widthInBlocks, heightInBlocks, depthMotion, numberOfBlocks;
	private ImageStack inStack, outStack;
	//private FloatProcessor outIp;
	private int blockSize = 8, blockSizeMotion = 32;
	
	int iterations = 5;
	
	// Current frame
	int curFrame = 2;
	
	// Vector alternatives for current frame
	private ArrayList<Vector3D> v_k;
	
	// Current vector field
	private Vector3D[] V_n;
	
	// Vector field of previous image
	private Vector3D[] V_n_previous;
	
	//
	private double g_l;
	

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_ALL | NO_CHANGES | STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {

		// dimensions stack
		width = ip.getWidth();
		widthInBlocks = width / blockSize;
		height = ip.getHeight();
		heightInBlocks = height / blockSize;
		inStack = imp.getStack();
		depth = inStack.getSize();
		depthMotion = depth - 1;
		numberOfBlocks = widthInBlocks * heightInBlocks;

		float[][] values = new float[depth][];
		Roi[] forms = new Roi[numberOfBlocks];
		Arrow tmp; // motion
		//OvalRoi tmp2; // no motion
		outStack = new ImageStack(width * 4, height * 4);
		FloatProcessor outIp;
		FloatProcessor Arrows;
		double xStart, yStart, xEnd, yEnd;
		int indexAlterMin;
		double [] costPerAlter = new double [v_k.size()];

		// read data from Stack
		for (int slice = 0; slice < depth; slice++) {
			values[slice] = (float[]) inStack.getProcessor(slice + 1).getPixels();
		}
//		System.out.println(v_k.size());
		//calculate and write motion vector field
		for (int slice = 0; slice < depthMotion; slice++) {
			for (int i = 0; i < iterations; i++) {
				for (int k = 0; k < numberOfBlocks; k++) {
					for (int a = 0; a < v_k.size(); a++) {
						costPerAlter [a] = costFunc(k, a);
					}
					indexAlterMin = minimizeCostFunc(k);
				}
			}
			
			Arrows = new FloatProcessor(width * 4, height * 4);
			outIp = new FloatProcessor(width * 4, height * 4);
			outIp.copyBits(inStack.getProcessor(slice+1).resize(4 * width), 0, 0, Blitter.COPY);
			for (int i = 0; i < heightInBlocks; i++) {
				for (int j = 0; j < widthInBlocks; j++) {
					
					//estimated motion vector
					//TO DO: +0.5 durch x und y koord des Bewegungsvektor ersetzen
					xStart = blockSizeMotion/ 2 + j * blockSizeMotion;
					yStart = blockSizeMotion / 2 + i * blockSizeMotion;
					xEnd = (blockSizeMotion/ 2 + j * blockSizeMotion) + 0.5;
					yEnd = (blockSizeMotion / 2 + i * blockSizeMotion) + 0.5;

					tmp = new Arrow(xEnd, yEnd, xStart, yStart);
					tmp.setStyle(Arrow.NOTCHED);
					forms[i] = tmp;

					// tmp2 = new OvalRoi(xStart, yStart, 6, 6);q
					// forms[i] = tmp2;
					Arrows.draw(forms[i]);
				}
			}
			Arrows.multiply(255/Arrows.maxValue());
			outIp.copyBits(Arrows, 0, 0, Blitter.COPY_ZERO_TRANSPARENT);
			outStack.addSlice(outIp);
		}

		String imgTitle = imp.getTitle();
		int l = imgTitle.length();
		int k = imgTitle.indexOf(" ");
		String title = k == -1 ? imgTitle : imgTitle.substring(60, l);
		ImagePlus window = new ImagePlus("Motion Vector Field of" + title, outStack);
		window.show();

	}
	
	/**
	 * Minimizes costFunc 
	 *
	 * @param  k	the index of the current block (starting with 0 in the upper left corner of the image).
	 * @return      the index of the vector with the lowest costs out of the 14 alternatives.
	 */
	private int minimizeCostFunc(int k){
		Iterator<Vector3D> alternatives = v_k.iterator();
		double min = Double.MAX_VALUE;
		double cur = Double.MAX_VALUE;
		int minIndex = -1;
		int index = 0;
		while(alternatives.hasNext()){
			cur = costFunc(k, index); 
			if(cur < min){
				min = cur;
				minIndex = index;
			}
			index++;
		}
		return minIndex;
	}
	
	private double costFunc(int k, int altIndex){
		
		return 0;
	}
	
	private double dataTerm (int k, int altIndex){
		
		return 0;
	}
	
	private double tempCoherence (int k, int altIndex){
		
		return 0;
	}
	
	private double spatialCoherence (int k, int altIndex){
		
		return 0;
	}
	
	private void getAlternatives (int k){
		
	}
}
