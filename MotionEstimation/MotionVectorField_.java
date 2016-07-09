import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Arrow;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class MotionVectorField_ implements PlugInFilter {

	private ImagePlus imp;
	private int width, height, depth, widthInBlocks, heightInBlocks, depthMotion, numberOfBlocks;
	private ImageStack inStack, outStack;
	private FloatProcessor outIp;
	private int blockSize = 8, blockSizeMotion=32;
	

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
		OvalRoi tmp2; // no motion
		outStack = new ImageStack(width*4, height*4);
		FloatProcessor outIp;
		double xStart, yStart, xEnd, yEnd;

		// read data from Stack
		for (int slice = 0; slice < depth; slice++) {
			values[slice] = (float[]) inStack.getProcessor(slice + 1).getPixels();
		}
		// write motion vector field
		for (int slice = 0; slice < depthMotion; slice++) {
			outIp = new FloatProcessor(width*4, height*4);

			for (int i = 0; i < heightInBlocks; i++) {
				for (int j = 0; j < widthInBlocks; j++) {
					
					xStart = blockSizeMotion/ 2 + j * blockSizeMotion;
					yStart = blockSizeMotion / 2 + i * blockSizeMotion;
					xEnd = (blockSizeMotion/ 2 + j * blockSizeMotion) + 0.5;
					yEnd = (blockSizeMotion / 2 + i * blockSizeMotion) + 0.5;

					tmp = new Arrow(xEnd, yEnd, xStart, yStart);
					tmp.setStyle("notched");
					forms[i] = tmp;

					// tmp2 = new OvalRoi(xStart, yStart, 6, 6);
					// forms[i] = tmp2;

					outIp.draw(forms[i]);
				}
			}
			outStack.addSlice(outIp);
		}

		String imgTitle = imp.getTitle();
		int l = imgTitle.length();
		int k = imgTitle.indexOf(" ");
		String title = k == -1 ? imgTitle : imgTitle.substring(60, l);
		ImagePlus window = new ImagePlus("Motion Vector Field of" + title, outStack);
		window.show();

	}
}
