import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

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
	// Inputs and output
	private ImagePlus imp;
	private ImageStack inStack, outStack;

	// Dimensions
	private int width, height, depth, widthInBlocks, heightInBlocks, depthMotion, numberOfBlocks,vecFieldWidth, vecFieldHeight;

	// Blocksizes
	private int blockSize = 8, blockSizeMotion = 32;

	// Numbers of Iterations
	private int iterations = 10;

	// Current frame
	private int curFrame = 0;

	// Exponent for temporal and spatial coherence
	private double nu = 1.3;

	// Lambda
	private double lambda;

	// Lambda_T
	private double lambda_T;

	// Picture data Y_n [frame][width * height]
	private float[][] values;

	// Vector alternatives for current frame
	private ArrayList<Vector2D> v_k;

	// Current vector field
	private Vector2D[] V_n;

	// Vector field of previous image
	private Vector2D[] V_n_previous;

	// Weighting factor
	private double g_l;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		lambda = 3.0;
		lambda_T = 3.0;
		v_k = new ArrayList<Vector2D>();
		return DOES_ALL | NO_CHANGES | STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {

		// Initializations
		width = ip.getWidth();
		widthInBlocks = width / blockSize;
		height = ip.getHeight();
		heightInBlocks = height / blockSize;
		inStack = imp.getStack();
		depth = inStack.getSize();
//		depth = 20;
		depthMotion = depth - 1;
		numberOfBlocks = widthInBlocks * heightInBlocks;
		vecFieldWidth = width * blockSizeMotion / blockSize;
		vecFieldHeight = height * blockSizeMotion / blockSize;

		values = new float[depth][];

		outStack = new ImageStack(vecFieldWidth , vecFieldHeight);
		FloatProcessor outIp;
		int indexAlterMin;

		// Read data from Stack
		for (int slice = 0; slice < depth; slice++) {
			values[slice] = (float[]) inStack.getProcessor(slice + 1).getPixels();
		}

		// Create a first vector field, where all vectors are equal (0, 0, 1)
		V_n_previous = new Vector2D[widthInBlocks * heightInBlocks];
		V_n = new Vector2D[widthInBlocks * heightInBlocks];
		for (int i = 0; i < widthInBlocks * heightInBlocks; i++) {
			V_n[i] = new Vector2D(0, 0);
		}

		// Calculate and write motion vector field
		for (int slice = 0; slice < depthMotion; slice++) {
			if(slice != 0){
				curFrame = slice;
				for (int i = 0; i < iterations; i++) {
					V_n_previous = V_n;
					for (int k = 0; k < numberOfBlocks; k++) {
						v_k = getAlternatives(k);
						indexAlterMin = minimizeCostFunc(k);
						V_n[k] = v_k.get(indexAlterMin);
					}
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
		Iterator<Vector2D> alternatives = v_k.iterator();
		double min = Double.MAX_VALUE;
		double cur = Double.MAX_VALUE;
		int minIndex = -1;
		int index = 0;

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
		double dataTerm = dataTerm(k, altIndex) * dataTerm(k, altIndex);
//		double dataTerm1 = Math.abs(dataTerm(k, altIndex));
		double spatialCoherence = lambda * spatialCoherence(k, altIndex);
		double tempCoherence = lambda_T * tempCoherence(k, altIndex);
//		System.out.println("dataTerm: " + dataTerm);
//		System.out.println("dataTerm1: " + dataTerm1);
//		System.out.println("Spatial: " + spatialCoherence);
//		System.out.println("Temporal: " + tempCoherence);

		return dataTerm + spatialCoherence + tempCoherence;
	}

	/**
	 * Calculates the data term of the cost function
	 *
	 * @param k   the index of the current block (starting with 0 in the upper
	 *            left corner of the image).
	 *
	 * @param altIndex
	 *            the index of the current alternative vector.
	 *
	 * @return the calculated value of the data term.
	 */
	 private double dataTerm(int k, int altIndex) {
 		// Current frame
 //		float[] Y_n = values[curFrame];
 		ImageProcessor Y_n = imp.getStack().getProcessor(curFrame + 1);

 		// Previous frame
 //		float[] Y_n_previous = values[curFrame - 1];

 		ImageProcessor Y_n_previous = imp.getStack().getProcessor(curFrame);

 		// Pixel starting positions calculated from k
 //		int posY = (k / widthInBlocks) * blockSize;
 //		int posX = k % widthInBlocks * blockSize;

 		double posY = (k / widthInBlocks) * blockSize;
 		double posX = k % widthInBlocks * blockSize;

 //		int posYprevious = (int) (posY + v_k.get(altIndex).getY());
 //		int posXprevious = (int) (posX + v_k.get(altIndex).getX());

 		double posYprevious = posY + v_k.get(altIndex).getY();
 		double posXprevious = posX + v_k.get(altIndex).getX();

 		// Pixel positions for cost calculation
 //		int posYcur = 0, posXcur = 0;
 //		int posYpreviousCur = 0, posXpreviousCur = 0;
 //		int imagePreviousAddress = 0;

 		double posYcur = 0, posXcur = 0;
 		double posYpreviousCur = 0, posXpreviousCur = 0;
 		double imagePreviousAddress = 0;

 		double result = 0;
 		for (int i = 0; i < blockSize * blockSize; i++) {
 			posYcur = posY + i / blockSize;
 			posXcur = posX + i % blockSize;
 			posYpreviousCur = posYprevious + i / blockSize;
 			posXpreviousCur = posXprevious + i % blockSize;

 			imagePreviousAddress = posYpreviousCur * width + posXpreviousCur;

 			// Edge protection (position + motion vector could overstep image boundaries)
 			if(imagePreviousAddress > width * height - 1)
 				imagePreviousAddress = width * height - 1;
 			if(imagePreviousAddress < 0)
 				imagePreviousAddress = 0;

 //			result += Y_n[posYcur * width + posXcur] - Y_n_previous[imagePreviousAddress];
 			result += Y_n.getInterpolatedValue(posXcur, posYcur) - Y_n_previous.getInterpolatedValue(posXpreviousCur, posYpreviousCur);
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

					Vector2D diff = v_k.get(altIndex).subtract(V_n[k + deltaX + deltaY * widthInBlocks]);

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
	private ArrayList<Vector2D> getAlternatives(int k) {

		ArrayList<Vector2D> nominees = new ArrayList<Vector2D>();

		// vector at the same position from the previous frame
		nominees.add(V_n_previous[k]);

		int[] indices = checkNeighbourIndices(k);

		double sumXcomponent = 0.0;
		double sumYcomponent = 0.0;
		double d = 0.0;

		for (int deltaX = indices[0]; deltaX <= indices[2]; deltaX++) {
			for (int deltaY = indices[1]; deltaY <= indices[3]; deltaY++) {

				if (!(deltaX == 0 && deltaY == 0)) {

					Vector2D neighbour = V_n[k + deltaX + deltaY * widthInBlocks];
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
		nominees.add(new Vector2D(sumXcomponent / d, sumYcomponent / d));

		// modification by 0.5 pixels
		nominees.add(new Vector2D(V_n[k].getX() + 0.5, V_n[k].getY()));
		nominees.add(new Vector2D(V_n[k].getX() - 0.5, V_n[k].getY()));
		nominees.add(new Vector2D(V_n[k].getX(), V_n[k].getY() + 0.5));
		nominees.add(new Vector2D(V_n[k].getX(), V_n[k].getY() - 0.5));

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
		int border=2;

		Arrows = new FloatProcessor(width * 4, height * 4);
		out = new FloatProcessor(width * 4, height * 4);
		out.copyBits(inStack.getProcessor(slice + 1).resize(4 * width), 0, 0, Blitter.COPY);
		for (int i = border; i < heightInBlocks-border; i++) {
			for (int j = border; j < widthInBlocks-border; j++) {

				x=V_n[i * (widthInBlocks) + j].getX();
				y=V_n[i * (widthInBlocks) + j].getY();

				// estimated motion vector
				xStart = blockSizeMotion / 2 + j * blockSizeMotion;
				yStart = blockSizeMotion / 2 + i * blockSizeMotion;
				xEnd = (blockSizeMotion / 2 + j * blockSizeMotion) + x;
				yEnd = (blockSizeMotion / 2 + i * blockSizeMotion) + y;

				if(Math.abs(x)+Math.abs(y)<=4.0){
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
