import java.util.ArrayList;
import java.util.Iterator;

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
	private int width, height, depth, widthInBlocks, heightInBlocks, depthMotion, numberOfBlocks,vecFieldWidth, vecFieldHeight;
	private ImageStack inStack, outStack;
	//Blocksizes
	private int blockSize = 8, blockSizeMotion = 32;

	// Numbers of Iterations
	private int iterations = 5;

	// Current frame
	private int curFrame = 0;

	// Exponent for temporal and spatial coherence
	private double nu = 1.3;

	// Lambda
	private double lambda;

	// Lambda_T
	private double lambda_T;

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
		lambda = 0.5;
		lambda_T = 0.5;
		v_k = new ArrayList<Vector3D>();
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
		vecFieldWidth=width * blockSizeMotion / blockSize;
		vecFieldHeight=height * blockSizeMotion / blockSize;

		float[][] values = new float[depth][];

		outStack = new ImageStack(vecFieldWidth , vecFieldHeight);
		FloatProcessor outIp;
		int indexAlterMin;

		// read data from Stack
		for (int slice = 0; slice < depth; slice++) {
			values[slice] = (float[]) inStack.getProcessor(slice + 1).getPixels();
		}

		// create a first vector field, where all vectors are equally (1,1,1)
		V_n_previous = new Vector3D[widthInBlocks * heightInBlocks];
		V_n = new Vector3D[widthInBlocks * heightInBlocks];
		for (int i = 0; i < widthInBlocks * heightInBlocks; i++) {
			V_n_previous[i] = new Vector3D(1, 1, 1);
			V_n[i] = new Vector3D(1, 1, 1);
		}
//		V_n = new Vector3D[widthInBlocks * heightInBlocks];
//		for (int i = 0; i < widthInBlocks * heightInBlocks; i++) {
//			V_n[0] = new Vector3D(1, 1, 1);
//		}

		// System.out.println(v_k.size());
		// calculate and write motion vector field
		for (int slice = 0; slice < depthMotion; slice++) {
			if(slice!=0){
				curFrame = slice;
				for (int i = 0; i < iterations; i++) {
					for (int k = 0; k < numberOfBlocks; k++) {
						v_k = getAlternatives(k);
						indexAlterMin = minimizeCostFunc(k);
						V_n[k] = v_k.get(indexAlterMin);
					}
					V_n_previous = V_n;
				}
			}

			outIp = generateMVFPerSlice(slice);
			outStack.addSlice(outIp);
		}

		String title = generateTitle();
		ImagePlus window = new ImagePlus("Motion Vector Field of " + title, outStack);
		window.show();

	}

	/**
	 * Minimizes costFunc
	 *
	 * @param k
	 *            the index of the current block (starting with 0 in the upper
	 *            left corner of the image).
	 * @return the index of the vector with the lowest costs out of the 14
	 *         alternatives.
	 */
	private int minimizeCostFunc(int k) {
		Iterator<Vector3D> alternatives = v_k.iterator();
		double min = Double.MAX_VALUE;
		double cur = Double.MAX_VALUE;
		int minIndex = -1;
		int index = 0;

//		for(Vector3D vec : v_k){				// use for each instead of iterator
		while (alternatives.hasNext()) {
			cur = costFunc(k, index);
			if (cur < min) {
				min = cur;
				minIndex = index;
			}
			alternatives.next();
			index++;
		}
		return minIndex;
	}

	/**
	 * Calculates the cost for one block
	 *
	 * @param k
	 *            the index of the current block (starting with 0 in the upper
	 *            left corner of the image).
	 * @param altIndex
	 *            the index of the current alternative vector.
	 * @return the calculated cost.
	 */
	private double costFunc(int k, int altIndex) {
		return 	dataTerm(k, altIndex) * dataTerm(k, altIndex)
				+ lambda * spatialCoherence(k, altIndex)
				+ lambda_T * tempCoherence(k, altIndex);
	}

	/**
	 * Calculates the data term of the cost function
	 *
	 * @param k
	 *            the index of the current block (starting with 0 in the upper
	 *            left corner of the image).
	 * @return the calculated value of the data term.
	 */
	private double dataTerm(int k, int altIndex) {
		float Y_n = values[curFrame][];
		float Y_n_previous = values[curFrame-1][];
		int posY = (k / widthInBlocks) * blockSize;
		int posX = k % widthInBlocks * blockSize;
		int posYprevious = (int) (posY + V_n[k].getY());
		int posXprevious = (int) (posX + V_n[k].getX());
		double result = 0;
		for (int i = 0; i < blockSize * blockSize; i++) {
			posY = posY + i / blockSize;
			posX = posX + i % blockSize;
			posYprevious = posYprevious + i / blockSize;
			posXprevious = posXprevious + i % blockSize;
			result += Y_n[posX * posY] - Y_n_previous[posXprevious * posYprevious];
		}
		return result;
	}

	/**
	 * Calculates the temporal coherence function
	 *
	 * @param k
	 *            the index of the current block (starting with 0 in the upper
	 *            left corner of the image).
	 * @param altIndex
	 *            the index of the current alternative vector.
	 * @return the calculated temporal coherence value.
	 */
	private double tempCoherence(int k, int altIndex) {
		return Math.pow(V_n_previous[k].distance(v_k.get(altIndex)), nu);
	}

	/**
	 * Calculates the spatial coherence/local energy
	 *
	 * @param k
	 *            the index of the current block
	 * @param altIndex
	 *            the index of the alternative vector
	 * @return value of calculated spatial coherence/local energy
	 */
	private double spatialCoherence(int k, int altIndex) {

		double res = 0;
		double n;

		int indices[] = checkNeighbourIndices(k);

		for (int deltaX = indices[0]; deltaX <= indices[2]; deltaX++) {
			for (int deltaY = indices[1]; deltaY <= indices[3]; deltaY++) {

				if (!(deltaX == 0 && deltaY == 0)) {

					Vector3D diff = v_k.get(altIndex).subtract(V_n[k + deltaX + deltaY * widthInBlocks]);

					n = Math.sqrt(diff.getX() * diff.getX() + diff.getY() * diff.getY());
					n = Math.pow(n, nu);

					if (deltaX == 0 || deltaY == 0)
						g_l = 1;
					else
						g_l = 0.5;

					res += g_l * n;
				}
			}
		}
		return res;
	}

	/**
	 * Determines the alternative vectors according to the paper by T. Aach and
	 * D. Kunz. Inner blocks will have 14 alternatives. Bordering blocks will
	 * have less alternatives depending on their location.
	 *
	 * @param k
	 *            the index of the current block
	 * @return An ArrayList containing the alternative vectors
	 */
	private ArrayList<Vector3D> getAlternatives(int k) {

		ArrayList<Vector3D> nominees = new ArrayList<Vector3D>();

		// vector at the same position from the previous frame
		nominees.add(V_n_previous[k]);

		int[] indices = checkNeighbourIndices(k);

		double sumXcomponent = 0.0;
		double sumYcomponent = 0.0;
		double d = 0.0;

		for (int deltaX = indices[0]; deltaX <= indices[2]; deltaX++) {
			for (int deltaY = indices[1]; deltaY <= indices[3]; deltaY++) {

				if (!(deltaX == 0 && deltaY == 0)) {

					Vector3D neighbour = V_n[k + deltaX + deltaY * widthInBlocks];
					nominees.add(neighbour);

					if (deltaX == 0 || deltaY == 0)
						g_l = 1;
					else
						g_l = 0.5;

					sumXcomponent += g_l * neighbour.getX();
					sumXcomponent += g_l * neighbour.getY();
					d += g_l;
				}
			}
		}

		// add vector of weighted average values
		nominees.add(new Vector3D(sumXcomponent / d, sumYcomponent / d, 1));

		// modification by 0.5 pixels
		nominees.add(new Vector3D(V_n[k].getX() + 0.5, V_n[k].getY(), 1));
		nominees.add(new Vector3D(V_n[k].getX() - 0.5, V_n[k].getY(), 1));
		nominees.add(new Vector3D(V_n[k].getX(), V_n[k].getY() + 0.5, 1));
		nominees.add(new Vector3D(V_n[k].getX(), V_n[k].getY() - 0.5, 1));

		return nominees;
	}

	/**
	 * Determines the indices in order to iterate over neighbouring blocks.
	 * Bordering blocks having less than eight neighbours will result in
	 * alternative indices.
	 *
	 * @param k
	 *            the index of the current block
	 * @return Array containing the indices
	 */
	private int[] checkNeighbourIndices(int k) {

		int[] neighbourIndices = new int[4];

		neighbourIndices[0] = -1; // startX
		neighbourIndices[1] = -1; // startY
		neighbourIndices[2] = 1; // endX
		neighbourIndices[3] = 1; // endY

		// first row
		if (k / widthInBlocks == 0) {
			neighbourIndices[0] = -1;
			neighbourIndices[1] = 0;
			neighbourIndices[2] = 1;
			neighbourIndices[3] = 1;

			// left upper corner
			if (k == 0)
				neighbourIndices[0] = 0;

			// right upper corner
			if (k == widthInBlocks - 1)
				neighbourIndices[2] = 0;
		}

		// last row
		if (k / widthInBlocks == heightInBlocks - 1) {
			neighbourIndices[0] = -1;
			neighbourIndices[1] = -1;
			neighbourIndices[2] = 1;
			neighbourIndices[3] = 0;

			// left lower corner
			if (k == 0)
				neighbourIndices[0] = 0;

			// right lower corner
			if (k == numberOfBlocks-1)				// 'numberOfBlocks' instead of 'widthInBlocks'
				neighbourIndices[2] = 0;
		}

		// left side
		if (k / widthInBlocks != 0 && k / widthInBlocks != heightInBlocks - 1 && k % widthInBlocks == 0) {
			neighbourIndices[0] = 0;
			neighbourIndices[1] = -1;
			neighbourIndices[2] = 1;
			neighbourIndices[3] = 1;
		}

		// right side
		if (k / widthInBlocks != 0 && k / widthInBlocks != heightInBlocks - 1
				&& k % widthInBlocks == widthInBlocks - 1) {
			neighbourIndices[0] = -1;
			neighbourIndices[1] = -1;
			neighbourIndices[2] = 0;
			neighbourIndices[3] = 1;
		}

		return neighbourIndices;
	}


	FloatProcessor generateMVFPerSlice(int slice){
		FloatProcessor Arrows,out;
		double xStart, yStart, xEnd, yEnd;
		Roi[] forms = new Roi[numberOfBlocks];
		Arrow tmp; // motion
		OvalRoi tmp2; // no motion
		double x,y;

		Arrows = new FloatProcessor(width * 4, height * 4);
		out = new FloatProcessor(width * 4, height * 4);
		out.copyBits(inStack.getProcessor(slice + 1).resize(4 * width), 0, 0, Blitter.COPY);
		for (int i = 0; i < heightInBlocks; i++) {
			for (int j = 0; j < widthInBlocks; j++) {

				x=V_n[i * widthInBlocks + j].getX();
				y=V_n[i * widthInBlocks + j].getY();

				// estimated motion vector
				xStart = blockSizeMotion / 2 + j * blockSizeMotion;
				yStart = blockSizeMotion / 2 + i * blockSizeMotion;
				xEnd = (blockSizeMotion / 2 + j * blockSizeMotion) + x;
				yEnd = (blockSizeMotion / 2 + i * blockSizeMotion) + y;

				if(Math.abs(x)+Math.abs(y)<=0.1){
					tmp2 = new OvalRoi(xEnd, yEnd, 6, 6);
					forms[i] = tmp2;
				}

				else{
					tmp = new Arrow(xEnd, yEnd, xStart, yStart);
					tmp.setStyle(Arrow.NOTCHED);
					forms[i] = tmp;
				}

				Arrows.draw(forms[i]);
			}
		}
		Arrows.multiply(255 / Arrows.maxValue());
		out.copyBits(Arrows, 0, 0, Blitter.COPY_ZERO_TRANSPARENT);
		return out;
	}

	String generateTitle(){
		String imgTitle = imp.getTitle();
		int l = imgTitle.length();
		int k = imgTitle.indexOf(" ");
		String t = k == -1 ? imgTitle : imgTitle.substring(60, l);
		return t;
	}
}
