/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;
import ij.measure.ResultsTable;
import ij.macro.Functions;
import ij.macro.MacroExtension;
import ij.macro.ExtensionDescriptor;

public class fit_RICS_jru_v1 implements PlugInFilter, NLLSfitinterface, MacroExtension
{
	ImagePlus imp;
	boolean c2test;
	int xc, yc, xpts, xskip, ypts, maxiter, photons, psfIndex, photonsIndex;
	int fullWidth, fullHeight;
	double pi = 3.14159265;
	double pixelsize, pixeltime, linetime;	
	int   [] fixes;
	float [] pixels;
	float [] g0val;
	float [] ac;
	float [] xaxis;
	float [] yaxis;
	float [][] acdisp;
	float [][] fitdisp;
	double[] params;
	double[] stats;
	double[][] constraints;
	String[] psfOptions     = { "3D Gaussian", "2D Gaussian" };
	String[] photonsOptions = { "mono-photon", "bi-photon" };
	Checkbox cb1;
	TextField tf1  = new TextField(10);

	public int setup(String arg, ImagePlus imp)
	{
		this.imp = imp;
		return DOES_ALL;
	}

	public void run(ImageProcessor ip)
	{
		pixels      = (float[]) ip.getPixels();
		fullWidth   = imp.getWidth();
		fullHeight  = imp.getHeight();
		
		if(!showDialog1())
			return;

		fixes       = new int[7];
		fixes[0]    = fixes[1] = 1;
		fixes[5]    = fixes[6] = 1;
		stats       = new double[2];
		constraints = new double[][] {{ 0.1, 0.0, -10.0, 0.0, 0.001, 0.0, 0.001 }, { 2.0, 20.0, 10.0, 1000.0, 1000.0, 1000.0, 1000.0 }};
		ac          = new float[xpts* ypts];
		acdisp      = new float[xpts][ypts];
		fitdisp     = new float[xpts][ypts];
		xaxis       = new float[xpts];
		yaxis       = new float[ypts];
		for(int i = 0; i < ypts; i++)
		{
			yaxis[i] = (float) pixelsize * (float) i;
			for(int j = 0; j < xpts; j++)
			{
				if(i == 0)
					xaxis[j] = (float) pixelsize * (float) j;
				ac[j + i * xpts] = pixels[(i + yc) * fullWidth + j + xc];
				acdisp[j][i]     = ac[j + i * xpts];
			}
		}
		g0val = new float[xskip];
		for(int i = 0;i < xskip; i++)
			g0val[i] = ac[i];

		if (IJ.macroRunning())
		{
			Functions.registerExtensions(this);
			fitting();
		}
		else
		{
			PlotWindow3D pw = new PlotWindow3D("RICS fit", "x (\u00B5m)", "y (\u00B5m)", "G(\u03BE,\u03C8)", xaxis, yaxis, acdisp);
			pw.draw();
			pw.addPoints(xaxis, yaxis, fitdisp, true);

			while(showDialog2())
			{
				fitting();
				pw.updateSeries(xaxis, yaxis, fitdisp, 1, true);
			}
		}
	}

	private boolean showDialog1()
	{
		//initialize the variables
		xc = fullWidth / 2; yc = fullHeight / 2; xpts = fullWidth / 2; ypts = fullHeight / 2; xskip = 2; maxiter = 10; c2test = false;
//		photons = 1; pixelsize = imp.getCalibration().pixelWidth; pixeltime = 0.0000028825; linetime = 0.00125;
		pixelsize = imp.getCalibration().pixelWidth; pixeltime = 0.00000288250; linetime = 0.0025;
				
		GenericDialog gd1 = new GenericDialog("Options");
		psfIndex     = 0;
		photonsIndex = 0;
		gd1.addChoice      ("PSF"             , psfOptions                     , psfOptions[psfIndex]);
		gd1.addNumericField("X_center"        , xc                             , 0);
		gd1.addNumericField("Y_center"        , yc                             , 0);
		gd1.addNumericField("Pts_in_x_to_skip", xskip                          , 0);
		gd1.addNumericField("Pts_in_x_to_fit" , xpts                           , 0);
		gd1.addNumericField("Pts_in_y_to_fit" , ypts                           , 0);
//		gd1.addNumericField("Photons"         , photons                        , 0);
		gd1.addChoice      ("Photons"         , photonsOptions                 , photonsOptions[0]);
		gd1.addNumericField("Pixel_Size"      , (float) (pixelsize)            , 7, 10, "\u00B5m");
		gd1.addNumericField("Pixel_Time"      , (float) (pixeltime * 1000000.0), 7, 10, "\u00B5s");
		gd1.addNumericField("Line_Time"       , (float) (linetime * 1000.0)    , 7, 10, "ms");
		gd1.showDialog();
		if(gd1.wasCanceled())
			return false;
		psfIndex     =       gd1.getNextChoiceIndex();
		xc           = (int) gd1.getNextNumber();
		yc           = (int) gd1.getNextNumber();
		xskip        = (int) gd1.getNextNumber();
		xpts         = (int) gd1.getNextNumber();
		ypts         = (int) gd1.getNextNumber();
//		photons      = (int) gd1.getNextNumber();
		photonsIndex =       gd1.getNextChoiceIndex();
		pixelsize    =      (gd1.getNextNumber());
		pixeltime    =      (gd1.getNextNumber()) / 1000000.0;
		linetime     =      (gd1.getNextNumber()) / 1000.0;
		photons      = photonsIndex + 1;
		params = new double[7];
		params[0] = 0.21; params[1] = 4.0;  params[2] = 1.0e-9;
		params[3] = 1.0 ; params[4] = 87.0; params[6] = 5.0;

		return true;
	}

	private boolean showDialog2()
	{
		GenericDialog gd2 = new GenericDialog("Fitting Parameters");
		gd2.setLayout(new GridLayout(12, 2));
		gd2.addCheckbox    ("Test_Chi_Squared"     , c2test);
		gd2.addMessage     ("");
		gd2.addMessage     ("");
		gd2.addNumericField("Max_Iterations"       , maxiter        ,  0, 15, " ");
		gd2.addMessage     ("");
		gd2.addNumericField("w0"                   , params[0]      , 10, 15, "\u00B5m");
		gd2.addCheckbox    ("fix_w0"               , (fixes[0] ==1 ));
		gd2.addNumericField("z0_/_w0"              , params[1]      , 10, 15, " ");
		gd2.addCheckbox    ("fix_z0_/_w0"          , (fixes[1] ==1 ));
		gd2.addNumericField("Baseline"             , params[2]      , 10, 15, " ");
		gd2.addCheckbox    ("fix_Baseline"         , (fixes[2] ==1 ));
		gd2.addNumericField("G(0)_1"               , params[3]      , 10, 15, " ");
		gd2.addCheckbox    ("fix_G(0)_1"           , (fixes[3] ==1 ));
		gd2.addNumericField("D1"                   , params[4]      , 10, 15, "\u00B5m^2 / sec");
		gd2.addCheckbox    ("fix_D1"               , (fixes[4] ==1 ));
		gd2.addNumericField("G(0)_2"               , params[5]      , 10, 15, " ");
		gd2.addCheckbox    ("fix_G(0)_2"           , (fixes[5] ==1 ));
		gd2.addNumericField("D2"                   , params[6]      , 10, 15, "\u00B5m^2 / sec");
		gd2.addCheckbox    ("fix_D2"               , (fixes[6] ==1 ));
		gd2.addNumericField("Iterations_Completed" , (int)stats[0]  , 10, 15, " ");
		gd2.addMessage     ("");		
		gd2.addNumericField("Chi_Squared"          , (float)stats[1], 10, 15, " ");
		gd2.addMessage     ("");
		gd2.addMessage     ("");
		gd2.showDialog();
		if(gd2.wasCanceled())
		{
			dataInResult();
//			dataInLog();
			return false;
		}
		c2test = gd2.getNextBoolean();
		maxiter = (int) gd2.getNextNumber();
		for(int i = 0; i < 7; i++)
		{
			params[i] = gd2.getNextNumber();
			fixes [i] = gd2.getNextBoolean() ? 1 : 0;
		}
		return true;
	}

	private void dataInResult()
	{
		ResultsTable results = new ResultsTable();
		results.incrementCounter();
		results.addValue ("Points in x to skip" , xskip);
		results.addValue ("Points in x to fit"  , xpts);
		results.addValue ("Points in y to fit"  , ypts);
		results.addValue ("Fix w0"              , intToBool(fixes[0]));
		results.addValue ("w0 (\u00B5m)"        , params[0]);
		results.addValue ("Fix z0/w0?"          , intToBool(fixes[1]));
		results.addValue ("z0/w0"               , params[1]);
		results.addValue ("Fix Baseline?"       , intToBool(fixes[2]));
		results.addValue ("Baseline"            , params[2]);
		results.addValue ("Fix G(0) 1?"         , intToBool(fixes[3]));
		results.addValue ("G(0) 1"              , params[3]);
		results.addValue ("Fix D1?"             , intToBool(fixes[4]));
		results.addValue ("D1 (\u00B5m^2 / sec)", params[4]);
		results.addValue ("Fix G(0) 2?"         , intToBool(fixes[5]));
		results.addValue ("G(0) 2"              , params[5]);
		results.addValue ("Fix D2?"             , intToBool(fixes[6]));
		results.addValue ("D2 (\u00B5m^2 / sec)", params[6]);
		results.show("Fit Results");
	}

	private void dataInLog()
	{
		IJ.log ("Fit Results");
		IJ.log ("Points in x to skip " + xskip);
		IJ.log ("Points in x to fit  " + xpts);
		IJ.log ("Points in y to fit  " + ypts);
		IJ.log ("Fix w0?             " + intToBool(fixes[0]));
		IJ.log ("w0 (\u00B5m)        " + (float) params[0]);
		IJ.log ("Fix z0/w0?          " + intToBool(fixes[1]));
		IJ.log ("z0/w0               " + (float) params[1]);
		IJ.log ("Fix Baseline?       " + intToBool(fixes[2]));
		IJ.log ("Baseline            " + (float) params[2]);
		IJ.log ("Fix G(0) 1?         " + intToBool(fixes[3]));
		IJ.log ("G(0) 1              " + (float) params[3]);
		IJ.log ("Fix D1?             " + intToBool(fixes[4]));
		IJ.log ("D1 (\u00B5m^2 / sec)" + (float) params[4]);
		IJ.log ("Fix G(0) 2?         " + intToBool(fixes[5]));
		IJ.log ("G(0) 2              " +(float) params[5]);
		IJ.log ("Fix D2?             " + intToBool(fixes[6]));
		IJ.log ("D2 (\u00B5m^2 / sec)" + (float) params[6]);
	}

	private String intToBool(int value)
	{
		return (value == 1) ? "true" : "false";
	}

	private void fitting()
	{
		NLLSfit nf;
		if(c2test)
			nf = new NLLSfit(this, 0.0001, 0, 0.1);				
		else
			nf = new NLLSfit(this, 0.0001, maxiter, 0.1);
		float[] fit = nf.fitdata(params, fixes, constraints, ac, null, stats, false);
		for(int i = 0;i < ypts; i++)
			for(int j = 0; j < xpts; j++)
				fitdisp[j][i] = fit[j + i * xpts];
	}

	public double fitfunc(double[] params, int indvar)
	{
		//the params list is w0 , z0 / w0, baseline, g01, D1, g02, D2
		if(indvar < xskip)
			return g0val[indvar];
		else
		{
			int    j            = (int)(indvar % xpts);
			int    i            = (int)((indvar - j) / xpts);
			double mult         = 4.0 * (double) photons;
			double xpixels      = (double) j;
			double xdistance    = xpixels * pixelsize;
			double ypixels      = (double) i;
			double ydistance    = ypixels * pixelsize;
			double tau          = pixeltime * xpixels + linetime * ypixels;
			double sqr_distance = xdistance * xdistance + ydistance * ydistance;
			double w02          = params[0] * params[0];
			double z02          = params[0] * params[1] * params[0] * params[1];
			double dumdouble1   = params[3] / (1.0 + ((mult * params[4] * tau) / w02));
			if(psfIndex == 0)
				dumdouble1 /= Math.sqrt(1.0 + ((mult * params[4] * tau) / z02));
			//note that the exponent differs from digman et al by a factor of 2
			dumdouble1 *= Math.exp(((double)photons) * (((-sqr_distance) / w02) / (1.0 + ((mult * params[4] * tau) / w02))));
			double dumdouble = dumdouble1;
			dumdouble1 = params[5] / (1.0 + ((mult * params[6] * tau) / w02));
			if(psfIndex == 0)
				dumdouble1 /= Math.sqrt(1.0 + ((mult * params[6] * tau) / z02));
			dumdouble1 *= Math.exp(((double) photons) * (((-sqr_distance) / w02) / (1.0 + ((mult * params[6] * tau) / w02))));
			return dumdouble1 + dumdouble + params[2];
		}
	}

	public void showresults(String results)
	{
		IJ.log(results);
	}

	private ExtensionDescriptor[] extensions =
	{
		ExtensionDescriptor.newDescriptor("getXpts"       , this),
		ExtensionDescriptor.newDescriptor("getXskip"      , this),
		ExtensionDescriptor.newDescriptor("getYpts"       , this),
		ExtensionDescriptor.newDescriptor("getFixw0"      , this),
		ExtensionDescriptor.newDescriptor("getW0"         , this),
		ExtensionDescriptor.newDescriptor("getFixz0_w0"   , this),
		ExtensionDescriptor.newDescriptor("getZ0_w0"      , this),
		ExtensionDescriptor.newDescriptor("getFixBaseline", this),
		ExtensionDescriptor.newDescriptor("getBaseline"   , this),
		ExtensionDescriptor.newDescriptor("getFixG0_1"    , this),
		ExtensionDescriptor.newDescriptor("getG0_1"       , this),
		ExtensionDescriptor.newDescriptor("getFixD1"      , this),
		ExtensionDescriptor.newDescriptor("getD1"         , this),
		ExtensionDescriptor.newDescriptor("getFixG0_2"    , this),
		ExtensionDescriptor.newDescriptor("getG0_2"       , this),
		ExtensionDescriptor.newDescriptor("getFixD2"      , this),
		ExtensionDescriptor.newDescriptor("getD2"         , this),
	};

	public ExtensionDescriptor[] getExtensionFunctions()
	{
		return extensions;
	}

	public String handleExtension(String name, Object[] args)
	{
		if (name.equals("getXpts"))
			return Integer.toString(xpts);
		else if (name.equals("getXskip"))
			return Integer.toString(xskip);
		else if (name.equals("getYpts"))
			return Integer.toString(ypts);
		else if (name.equals("getFixw0"))
			return intToBool(fixes[0]);
		else if (name.equals("getWw0"))
			return Double.toString(params[0]);
		else if (name.equals("getFixz0_w0"))
			return intToBool(fixes[1]);
		else if (name.equals("getZ0_w0"))
			return Double.toString(params[1]);
		else if (name.equals("getFixBaseline"))
			return intToBool(fixes[2]);
		else if (name.equals("getBaseline"))
			return Double.toString(params[2]);
		else if (name.equals("getFixG0_1"))
			return intToBool(fixes[3]);
		else if (name.equals("getG0_1"))
			return Double.toString(params[3]);
		else if (name.equals("getFixD1"))
			return intToBool(fixes[4]);
		else if (name.equals("getD1"))
			return Double.toString(params[4]);
		else if (name.equals("getFixG0_2"))
			return intToBool(fixes[5]);
		else if (name.equals("getG0_2"))
			return Double.toString(params[5]);
		else if (name.equals("getFixD2"))
			return intToBool(fixes[6]);
		else if (name.equals("getD2"))
			return Double.toString(params[6]);

		return null;
	}
}
