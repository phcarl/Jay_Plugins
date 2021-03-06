/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.util.Arrays;

public class algutils{

	/*********************************
	 * this returns my own array type index: 0 for byte, 1 for short, 2 for float, 3 for double, and 4 for int
	 * @param arr
	 * @return
	 */
	public static int get_array_type(Object arr){
		if(arr instanceof byte[])
			return 0;
		if(arr instanceof short[])
			return 1;
		if(arr instanceof float[])
			return 2;
		if(arr instanceof double[])
			return 3;
		if(arr instanceof int[])
			return 4;
		return -1;
	}
	
	/*********************************
	 * this gets the length of an array of undefined type
	 * @param arr
	 * @return
	 */
	public static int get_array_length(Object arr){
		int atype=get_array_type(arr);
		if(atype==0) return ((byte[])arr).length;
		if(atype==1) return ((short[])arr).length;
		if(atype==2) return ((float[])arr).length;
		if(atype==3) return ((double[])arr).length;
		if(atype==4) return ((int[])arr).length;
		return -1;
	}

	/************************************
	 * this converts an array to a specific type, making a copy.  See above for types
	 * @param oldarr
	 * @param type
	 * @return
	 */
	public static Object convert_array(Object oldarr,int type){
		// can convert from and to byte, short, or float
		switch(type){
		case 0:
			return convert_arr_byte(oldarr);
		case 1:
			return convert_arr_short(oldarr);
		case 2:
			return convert_arr_float(oldarr);
		case 3:
			return convert_arr_double(oldarr);
		case 4:
			return convert_arr_int(oldarr);
		}
		return null;
	}

	/*************************************
	 * makes a copy of an array as a byte
	 * @param oldarr
	 * @return
	 */
	public static byte[] convert_arr_byte(Object oldarr){
		if(oldarr instanceof double[]){
			double[] temparr=(double[])oldarr;
			byte[] newarr=new byte[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=(int)temparr[i];
				if(temp>=0){
					if(temp<256){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)255;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			byte[] newarr=new byte[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=(int)temparr[i];
				if(temp>=0){
					if(temp<256){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)255;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			byte[] newarr=new byte[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xffff;
				if(temp>=0){
					if(temp<256){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)255;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			byte[] newarr=new byte[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		if(oldarr instanceof int[]){
			int[] temparr=(int[])oldarr;
			byte[] newarr=new byte[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i];
				if(temp>=0){
					if(temp<256){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)255;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		return null;
	}
	
	public static byte[][] convert_arr_byte(Object[] input){
		byte[][] out=new byte[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=convert_arr_byte(input[i]);
		}
		return out;
	}

	/*****************************
	 * makes a copy of an array as a float
	 * @param oldarr
	 * @return
	 */
	public static float[] convert_arr_float(Object oldarr){
		if(oldarr instanceof double[]){
			double[] temparr=(double[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=(float)temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			float[] newarr=new float[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xffff;
				newarr[i]=temp;
			}
			return newarr;
		}
		if(oldarr instanceof int[]){
			int[] temparr=(int[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof long[]){
			long[] temparr=(long[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xff;
				newarr[i]=(float)temp;
			}
			return newarr;
		}
		return null;
	}
	
	public static float[][] convert_arr_float(Object[] input){
		float[][] out=new float[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=convert_arr_float(input[i]);
		}
		return out;
	}
	
	public static float[][] convert_arr_float2(Object[] input){
		float[][] out=new float[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=convert_arr_float2(input[i]);
		}
		return out;
	}

	/*****************************
	 * makes a copy of an array as an int
	 * @param oldarr
	 * @return
	 */
	public static int[] convert_arr_int(Object oldarr){
		if(oldarr instanceof double[]){
			double[] temparr=(double[])oldarr;
			int[] newarr=new int[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=(int)temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			int[] newarr=new int[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=(int)temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			int[] newarr=new int[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=temparr[i]&0xffff;
			}
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			int[] newarr=new int[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=temparr[i]&0xff;
			}
			return newarr;
		}
		if(oldarr instanceof int[]){
			int[] temparr=(int[])oldarr;
			int[] newarr=new int[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		return null;
	}
	
	public static int[][] convert_arr_int(Object[] input){
		int[][] out=new int[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=convert_arr_int(input[i]);
		}
		return out;
	}

	/*****************************
	 * makes a copy of an array as a double
	 * @param oldarr
	 * @return
	 */
	public static double[] convert_arr_double(Object oldarr){
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			double[] newarr=new double[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof double[]){
			double[] temparr=(double[])oldarr;
			double[] newarr=new double[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			double[] newarr=new double[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xffff;
				newarr[i]=temp;
			}
			return newarr;
		}
		if(oldarr instanceof int[]){
			int[] temparr=(int[])oldarr;
			double[] newarr=new double[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			double[] newarr=new double[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xff;
				newarr[i]=temp;
			}
			return newarr;
		}
		return null;
	}
	
	public static double[][] convert_arr_double(Object[] input){
		double[][] out=new double[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=convert_arr_double(input[i]);
		}
		return out;
	}

	/*****************************
	 * makes a copy of an array as a short
	 * @param oldarr
	 * @return
	 */
	public static short[] convert_arr_short(Object oldarr){
		if(oldarr instanceof double[]){
			double[] temparr=(double[])oldarr;
			short[] newarr=new short[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=(int)temparr[i];
				if(temp>=0){
					if(temp<65536){
						newarr[i]=(short)temp;
					}else{
						newarr[i]=(short)65536;
					}
				}else{
					newarr[i]=(short)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			short[] newarr=new short[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=(int)temparr[i];
				if(temp>=0){
					if(temp<65536){
						newarr[i]=(short)temp;
					}else{
						newarr[i]=(short)65536;
					}
				}else{
					newarr[i]=(short)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof int[]){
			int[] temparr=(int[])oldarr;
			short[] newarr=new short[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i];
				if(temp>=0){
					if(temp<65536){
						newarr[i]=(short)temp;
					}else{
						newarr[i]=(short)65536;
					}
				}else{
					newarr[i]=(short)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			short[] newarr=new short[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			short[] newarr=new short[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xff;
				newarr[i]=(short)temp;
			}
			return newarr;
		}
		return null;
	}
	
	public static short[][] convert_arr_short(Object[] input){
		short[][] out=new short[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=convert_arr_short(input[i]);
		}
		return out;
	}

	/*****************************
	 * converts an array to byte, doesn't copy if already a byte
	 * @param oldarr
	 * @return
	 */
	public static byte[] convert_arr_byte2(Object arr){
		if(arr instanceof byte[])
			return (byte[])arr;
		else
			return convert_arr_byte(arr);
	}

	/*****************************
	 * converts an array to short, doesn't copy if already a short
	 * @param oldarr
	 * @return
	 */
	public static short[] convert_arr_short2(Object arr){
		if(arr instanceof short[])
			return (short[])arr;
		else
			return convert_arr_short(arr);
	}

	/*****************************
	 * converts an array to float, doesn't copy if already a float
	 * @param oldarr
	 * @return
	 */
	public static float[] convert_arr_float2(Object arr){
		if(arr instanceof float[])
			return (float[])arr;
		else
			return convert_arr_float(arr);
	}

	/*****************************
	 * converts an array to int, doesn't copy if already a int
	 * @param oldarr
	 * @return
	 */
	public static int[] convert_arr_int2(Object arr){
		if(arr instanceof int[])
			return (int[])arr;
		else
			return convert_arr_int(arr);
	}

	/*****************************
	 * converts an array to double, doesn't copy if already a double
	 * @param oldarr
	 * @return
	 */
	public static double[] convert_arr_double2(Object arr){
		if(arr instanceof double[])
			return (double[])arr;
		else
			return convert_arr_double(arr);
	}
	
	public static float getPixelVal(Object arr,int x,int y,int width,int height,int dtype){
		if(x<0) return 0.0f;
		if(y<0) return 0.0f;
		if(x>=width) return 0.0f;
		if(y>=height) return 0.0f;
		return getArrVal(arr,x+y*width,dtype);
	}
	
	public static float getArrVal(Object arr,int index,int dtype){
		if(dtype==0) return ((byte[])arr)[index]&0xff;
		if(dtype==1) return ((short[])arr)[index]&0xffff;
		if(dtype==2) return ((float[])arr)[index];
		if(dtype==3) return (float)((double[])arr)[index];
		if(dtype==4) return ((int[])arr)[index];
		return Float.NaN;
	}

	/**********************
	 * converts a 1D image array to and array of rows
	 * @param arr: the array
	 * @param width
	 * @param height
	 * @return a 2D array with [height][width] indices
	 */
	public static float[][] array2multidim(float[] arr,int width,int height){
		float[][] temp=new float[height][width];
		int temp2=0;
		for(int i=0;i<height;i++){
			System.arraycopy(arr,temp2,temp[i],0,width);
			temp2+=width;
		}
		return temp;
	}

	/*******************
	 * copies a 2D float array
	 * @param arr
	 * @return
	 */
	public static float[][] clone_multidim_array(float[][] arr){
		float[][] temp=new float[arr.length][];
		for(int i=0;i<arr.length;i++){
			temp[i]=arr[i].clone();
		}
		return temp;
	}
	
	/*******************
	 * copies a 2D int array
	 * @param arr
	 * @return
	 */
	public static int[][] clone_multidim_array(int[][] arr){
		int[][] temp=new int[arr.length][];
		for(int i=0;i<arr.length;i++){
			temp[i]=arr[i].clone();
		}
		return temp;
	}

	/************************
	 * copies a 2D short array
	 * @param arr
	 * @return
	 */
	public static short[][] clone_multidim_array(short[][] arr){
		short[][] temp=new short[arr.length][];
		for(int i=0;i<arr.length;i++){
			temp[i]=arr[i].clone();
		}
		return temp;
	}

	/************************
	 * copies a 2D byte array
	 * @param arr
	 * @return
	 */
	public static byte[][] clone_multidim_array(byte[][] arr){
		byte[][] temp=new byte[arr.length][];
		for(int i=0;i<arr.length;i++){
			temp[i]=arr[i].clone();
		}
		return temp;
	}

	public static float getNeighborsAvg(Object image,int x,int y,int width,int height){
		Object neighbors=getNeighbors(image,x,y,width,height);
		if(neighbors!=null){
			return jstatistics.getstatistic("Avg",neighbors,null);
		}else{
			return 0.0f;
		}
	}

	public static float getNeighbors2Avg(Object image,int x,int y,int width,int height){
		Object neighbors=getNeighbors2(image,x,y,width,height);
		if(neighbors!=null){
			return jstatistics.getstatistic("Avg",neighbors,null);
		}else{
			return 0.0f;
		}
	}

	public static float getNeighbors2Stat(Object image,int x,int y,int width,int height,String stat){
		Object neighbors=getNeighbors2(image,x,y,width,height);
		if(neighbors!=null){
			return jstatistics.getstatistic(stat,neighbors,null);
		}else{
			return 0.0f;
		}
	}

	public static float getNeighborsAvg(Object[] image,int x,int y,int z,int width,int height){
		Object neighbors=getNeighbors(image,x,y,z,width,height);
		if(neighbors!=null){
			return jstatistics.getstatistic("Avg",neighbors,null);
		}else{
			return 0.0f;
		}
	}

	public static float getNeighbors2Avg(Object[] image,int x,int y,int z,int width,int height){
		Object neighbors=getNeighbors2(image,x,y,z,width,height);
		if(neighbors!=null){
			return jstatistics.getstatistic("Avg",neighbors,null);
		}else{
			return 0.0f;
		}
	}
	
	public static float getNeighbors2Stat(Object[] image,int x,int y,int z,int width,int height,String stat){
		Object neighbors=getNeighbors2(image,x,y,z,width,height);
		if(neighbors!=null){
			return jstatistics.getstatistic(stat,neighbors,null);
		}else{
			return 0.0f;
		}
	}

	public static Object getNeighbors(Object image,int x,int y,int width,int height){
		if(image instanceof float[]){
			return getNeighbors((float[])image,x,y,width,height);
		}else{
			if(image instanceof short[]) return getNeighbors((short[])image,x,y,width,height);
			else return getNeighbors((int[])image,x,y,width,height);
		}
	}

	public static Object getNeighbors2(Object image,int x,int y,int width,int height){
		if(image instanceof float[]){
			return getNeighbors2((float[])image,x,y,width,height);
		}else{
			return getNeighbors2((short[])image,x,y,width,height);
		}
	}

	public static Object getNeighbors(Object[] image,int x,int y,int z,int width,int height){
		if(z<=0||z>=(image.length-1)){
			return null;
		}
		if(x<=0||x>=(width-1)){
			return null;
		}
		if(y<=0||y>=(height-1)){
			return null;
		}
		boolean isfloat=(image[0] instanceof float[]);
		Object temp0=getNeighbors2(image[z-1],x,y,width,height);
		Object temp1=getNeighbors(image[z],x,y,width,height);
		Object temp2=getNeighbors2(image[z+1],x,y,width,height);
		if(isfloat){
			float[] temp=new float[9+9+8];
			System.arraycopy(temp0,0,temp,0,9);
			System.arraycopy(temp1,0,temp,9,8);
			System.arraycopy(temp2,0,temp,17,9);
			return temp;
		}else{
			short[] temp=new short[9+9+8];
			System.arraycopy(temp0,0,temp,0,9);
			System.arraycopy(temp1,0,temp,9,8);
			System.arraycopy(temp2,0,temp,17,9);
			return temp;
		}
	}

	public static Object getNeighbors2(Object[] image,int x,int y,int z,int width,int height){
		if(z<0||z>(image.length-1)){
			return null;
		}
		if(x<0||x>(width-1)){
			return null;
		}
		if(y<0||y>(height-1)){
			return null;
		}
		boolean isfloat=(image[0] instanceof float[]);
		Object temp0=null;
		Object temp1=null;
		Object temp2=null;
		int size=9;
		if(z>0){
			temp0=getNeighbors2(image[z-1],x,y,width,height);
			size+=9;
		}
		temp1=getNeighbors2(image[z],x,y,width,height);
		if(z<image.length){
			temp2=getNeighbors2(image[z+1],x,y,width,height);
			size+=9;
		}
		if(isfloat){
			float[] temp=new float[size];
			int counter=0;
			if(z>0){
				System.arraycopy(temp0,0,temp,0,9);
				counter+=9;
			}
			System.arraycopy(temp1,0,temp,counter,9);
			counter+=9;
			if(z<image.length)
				System.arraycopy(temp2,0,temp,counter,9);
			return temp;
		}else{
			short[] temp=new short[size];
			int counter=0;
			if(z>0){
				System.arraycopy(temp0,0,temp,0,9);
				counter+=9;
			}
			System.arraycopy(temp1,0,temp,counter,9);
			counter+=9;
			if(z<image.length)
				System.arraycopy(temp2,0,temp,counter,9);
			return temp;
		}
	}

	public static float[] getNeighbors(float[] objects,int x,int y,int width,int height){
		// here we don't return the center pixel and edges are returned null
		if(x<=0||x>=(width-1)){
			return null;
		}
		if(y<=0||y>=(height-1)){
			return null;
		}
		float[] temp=new float[8];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2+=2;
		temp[4]=objects[temp2];
		temp2+=(width-2);
		temp[5]=objects[temp2];
		temp2++;
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		return temp;
	}

	public static float[] getNeighbors2(float[] objects,int x,int y,int width,int height){
		// here we include the center pixel and edges are returned as possible
		if(x<0||x>(width-1)){
			return null;
		}
		if(y<0||y>(height-1)){
			return null;
		}
		if(x==0){
			if(y==0){
				// upper left
				float[] temp={objects[0],objects[1],objects[width],objects[width+1]};
				return temp;
			}else{
				if(y==(height-1)){
					// lower left
					int index=width*(height-1);
					float[] temp={objects[index],objects[index+1],objects[index-width],objects[index-width+1]};
					return temp;
				}else{
					// left side
					int index=(y-1)*width;
					float[] temp={objects[index],objects[index+1],objects[index+width],objects[index+width+1],objects[index+width+width],objects[index+width+width+1]};
					return temp;
				}
			}
		}
		if(x==(width-1)){
			if(y==0){
				// upper right
				float[] temp={objects[width-2],objects[width-1],objects[width+width-2],objects[width+width-1]};
				return temp;
			}else{
				if(y==(height-1)){
					// lower right
					int index=width*height-1;
					float[] temp={objects[index],objects[index-1],objects[index-width],objects[index-width-1]};
					return temp;
				}else{
					// right side
					int index=y*width-1;
					float[] temp={objects[index],objects[index-1],objects[index+width],objects[index+width-1],objects[index+width+width],objects[index+width+width-1]};
					return temp;
				}
			}
		}
		if(y==0){
			// already returned if in a corner
			// must be top edge
			int index=x-1;
			float[] temp={objects[index],objects[index+1],objects[index+2],objects[index+width],objects[index+width+1],objects[index+width+2]};
			return temp;
		}
		if(y==(height-1)){
			// must be on bottom edge
			int index=x-1+width*(height-1);
			float[] temp={objects[index],objects[index+1],objects[index+2],objects[index-width],objects[index-width+1],objects[index-width+2]};
			return temp;
		}
		float[] temp=new float[9];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2++;
		temp[4]=objects[temp2];
		temp2++;
		temp[5]=objects[temp2];
		temp2+=(width-2);
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		temp2++;
		temp[8]=objects[temp2];
		return temp;
	}

	public static short[] getNeighbors(short[] objects,int x,int y,int width,int height){
		if(x<=0||x>=(width-1)){
			return null;
		}
		if(y<=0||y>=(height-1)){
			return null;
		}
		short[] temp=new short[8];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2+=2;
		temp[4]=objects[temp2];
		temp2+=(width-2);
		temp[5]=objects[temp2];
		temp2++;
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		return temp;
	}
	
	public static int[] getNeighbors(int[] objects,int x,int y,int width,int height){
		if(x<=0||x>=(width-1)){
			return null;
		}
		if(y<=0||y>=(height-1)){
			return null;
		}
		int[] temp=new int[8];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2+=2;
		temp[4]=objects[temp2];
		temp2+=(width-2);
		temp[5]=objects[temp2];
		temp2++;
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		return temp;
	}

	public static short[] getNeighbors2(short[] objects,int x,int y,int width,int height){
		if(x<0||x>(width-1)){
			return null;
		}
		if(y<0||y>(height-1)){
			return null;
		}
		if(x==0){
			if(y==0){
				// upper left
				short[] temp={objects[0],objects[1],objects[width],objects[width+1]};
				return temp;
			}else{
				if(y==(height-1)){
					// lower left
					int index=width*(height-1);
					short[] temp={objects[index],objects[index+1],objects[index-width],objects[index-width+1]};
					return temp;
				}else{
					// left side
					int index=(y-1)*width;
					short[] temp={objects[index],objects[index+1],objects[index+width],objects[index+width+1],objects[index+width+width],objects[index+width+width+1]};
					return temp;
				}
			}
		}
		if(x==(width-1)){
			if(y==0){
				// upper right
				short[] temp={objects[width-2],objects[width-1],objects[width+width-2],objects[width+width-1]};
				return temp;
			}else{
				if(y==(height-1)){
					// lower right
					int index=width*height-1;
					short[] temp={objects[index],objects[index-1],objects[index-width],objects[index-width-1]};
					return temp;
				}else{
					// right side
					int index=y*width-1;
					short[] temp={objects[index],objects[index-1],objects[index+width],objects[index+width-1],objects[index+width+width],objects[index+width+width-1]};
					return temp;
				}
			}
		}
		if(y==0){
			// already returned if in a corner
			// must be top edge
			int index=x-1;
			short[] temp={objects[index],objects[index+1],objects[index+2],objects[index+width],objects[index+width+1],objects[index+width+2]};
			return temp;
		}
		if(y==(height-1)){
			// must be on bottom edge
			int index=x-1+width*(height-1);
			short[] temp={objects[index],objects[index+1],objects[index+2],objects[index-width],objects[index-width+1],objects[index-width+2]};
			return temp;
		}
		short[] temp=new short[9];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2++;
		temp[4]=objects[temp2];
		temp2++;
		temp[5]=objects[temp2];
		temp2+=(width-2);
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		temp2++;
		temp[8]=objects[temp2];
		return temp;
	}

	public static float get_region_stat(String stat,Object image,int x,int y,int rwidth,int rheight,int width,int height){
		float[] reg=get_region(image,x,y,rwidth,rheight,width,height);
		return jstatistics.getstatistic(stat,reg,null);
	}

	public static float[] get_region_stat(String stat,Object[] stack,int xc,int yc,int rwidth,int rheight,int width,int height){
		float[][] reg=get_region(stack,xc,yc,rwidth,rheight,width,height);
		float[] stats=new float[reg.length];
		for(int i=0;i<stats.length;i++)
			stats[i]=jstatistics.getstatistic(stat,reg[i],null);
		return stats;
	}

	/*****************
	 * gets a region of a stack
	 * @param stack
	 * @param xc: x center of the region
	 * @param yc: y center of the region
	 * @param rwidth
	 * @param rheight
	 * @param width
	 * @param height
	 * @return
	 */
	public static float[][] get_region(Object[] stack,int xc,int yc,int rwidth,int rheight,int width,int height){
		float[][] profile=new float[stack.length][];
		for(int i=0;i<stack.length;i++){
			profile[i]=get_region(stack[i],xc,yc,rwidth,rheight,width,height);
		}
		return profile;
	}
	
	/*****************
	 * gets a region of a stack
	 * @param stack
	 * @param x: x origin of the region
	 * @param y: y origin of the region
	 * @param rwidth
	 * @param rheight
	 * @param width
	 * @param height
	 * @return
	 */
	public static float[][] get_region2(Object[] stack,int x,int y,int rwidth,int rheight,int width,int height){
		float[][] profile=new float[stack.length][];
		for(int i=0;i<stack.length;i++){
			profile[i]=get_region2(stack[i],x,y,rwidth,rheight,width,height);
		}
		return profile;
	}

	public static float[] get_region(Object image,int xc,int yc,int rwidth,int rheight,int width,int height){
		// note that x and y are centers, not the origin
		if(image instanceof float[]){
			return get_region((float[])image,xc,yc,rwidth,rheight,width,height);
		}else{
			if(image instanceof byte[]){
				return get_region((byte[])image,xc,yc,rwidth,rheight,width,height);
			}else{
				return get_region((short[])image,xc,yc,rwidth,rheight,width,height);
			}
		}
	}
	
	public static float[] get_region2(Object image,int x,int y,int rwidth,int rheight,int width,int height){
		// note that x and y are the origin
		if(image instanceof float[]){
			return get_region2((float[])image,x,y,rwidth,rheight,width,height);
		}else{
			if(image instanceof byte[]){
				return get_region2((byte[])image,x,y,rwidth,rheight,width,height);
			}else{
				return get_region2((short[])image,x,y,rwidth,rheight,width,height);
			}
		}
	}

	public static float[] get_region(float[] image,int xc,int yc,int rwidth,int rheight,int width,int height){
		// here we return the rectangular region centered on x and y
		if(!inbounds(xc,yc,width,height)) return null;
		int startx=xc-rwidth/2;
		int starty=yc-rheight/2;
		return get_region2(image,startx,starty,rwidth,rheight,width,height);
	}

	public static float[] get_region(short[] image,int xc,int yc,int rwidth,int rheight,int width,int height){
		// here we return the rectangular region centered on x and y
		if(!inbounds(xc,yc,width,height)) return null;
		int startx=xc-rwidth/2;
		int starty=yc-rheight/2;
		return get_region2(image,startx,starty,rwidth,rheight,width,height);
	}

	public static float[] get_region(byte[] image,int xc,int yc,int rwidth,int rheight,int width,int height){
		// here we return the rectangular region centered on x and y
		if(!inbounds(xc,yc,width,height)) return null;
		int startx=xc-rwidth/2;
		int starty=yc-rheight/2;
		return get_region2(image,startx,starty,rwidth,rheight,width,height);
	}
	
	public static float[] get_region2(float[] image,int x,int y,int rwidth,int rheight,int width,int height){
		//here we return the region starting at x,y
		if(!inbounds(x,y,width,height)) return null;
		int startx=x; int starty=y;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		int xpix=endx-startx+1;
		int ypix=endy-starty+1;
		float[] output=new float[xpix*ypix];
		int counter=0;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				output[counter]=image[j+i*width];
				counter++;
			}
		}
		return output;
	}
	
	public static float[] get_region2(short[] image,int x,int y,int rwidth,int rheight,int width,int height){
		//here we return the region starting at x,y
		if(!inbounds(x,y,width,height)) return null;
		int startx=x; int starty=y;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		int xpix=endx-startx+1;
		int ypix=endy-starty+1;
		float[] output=new float[xpix*ypix];
		int counter=0;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				output[counter]=image[j+i*width]&0xffff;
				counter++;
			}
		}
		return output;
	}
	
	public static float[] get_region2(byte[] image,int x,int y,int rwidth,int rheight,int width,int height){
		//here we return the region starting at x,y
		if(!inbounds(x,y,width,height)) return null;
		int startx=x; int starty=y;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		int xpix=endx-startx+1;
		int ypix=endy-starty+1;
		float[] output=new float[xpix*ypix];
		int counter=0;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				output[counter]=image[j+i*width]&0xff;
				counter++;
			}
		}
		return output;
	}

	public static float[] get_region_pad(float[] image,int x,int y,int rwidth,int rheight,int width,int height){
		// here we return the rectangular region centered on x and y
		if(!inbounds(x,y,width,height))
			return null;
		int startx=x-rwidth/2;
		int starty=y-rheight/2;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		float[] output=new float[rwidth*rheight];
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				output[j-startx+(i-starty)+rwidth]=image[j+i*width];
			}
		}
		return output;
	}

	public static void set_region(float[] image,int x,int y,int rwidth,int rheight,int width,int height,float value){
		// here we set the rectangular region centered on x and y equal to value
		if(!inbounds(x,y,width,height))
			return;
		int startx=x-rwidth/2;
		int starty=y-rheight/2;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				image[j+i*width]=value;
			}
		}
	}
	
	public static void set_region(float[] image,int x,int y,int rwidth,int rheight,int width,int height,float[] value){
		// here we set the rectangular region centered on x and y equal to value
		if(!inbounds(x,y,width,height))
			return;
		int startx=x-rwidth/2;
		int starty=y-rheight/2;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		int counter=0;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				image[j+i*width]=value[counter];
				counter++;
			}
		}
	}

	public static void set_region(short[] image,int x,int y,int rwidth,int rheight,int width,int height,short value){
		// here we set the rectangular region centered on x and y equal to value
		if(!inbounds(x,y,width,height))
			return;
		int startx=x-rwidth/2;
		int starty=y-rheight/2;
		int endx=startx+rwidth-1;
		int endy=starty+rheight-1;
		if(startx<0)
			startx=0;
		if(starty<0)
			starty=0;
		if(endx>(width-1))
			endx=width-1;
		if(endy>(height-1))
			endy=height-1;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				image[j+i*width]=value;
			}
		}
	}

	public static boolean inbounds(int x,int y,int width,int height){
		if(x<0)
			return false;
		if(y<0)
			return false;
		if(x>(width-1))
			return false;
		if(y>(height-1))
			return false;
		return true;
	}

	public static void copy_subarray(Object source,int soff,Object dest,int doff,int length){
		//assume that source and dest have same data type
		if(source instanceof float[]){
			System.arraycopy(source,soff,dest,doff,length);
		}else{
			if(source instanceof short[]){
				System.arraycopy(source,soff,dest,doff,length);
			}else{
				if(source instanceof byte[]){
					System.arraycopy(source,soff,dest,doff,length);
				}else{
					System.arraycopy(source,soff,dest,doff,length);
				}
			}
		}
	}

	public static Object get_subarray(Object source,int off,int length){
		if(source instanceof float[]){
			float[] temp=new float[length];
			System.arraycopy(source,off,temp,0,length);
			return temp;
		}else{
			if(source instanceof short[]){
				short[] temp=new short[length];
				System.arraycopy(source,off,temp,0,length);
				return temp;
			}else{
				if(source instanceof byte[]){
					byte[] temp=new byte[length];
					System.arraycopy(source,off,temp,0,length);
					return temp;
				}else{
					int[] temp=new int[length];
					System.arraycopy(source,off,temp,0,length);
					return temp;
				}
			}
		}
	}

	public static Object clone_obj_array(Object source){
		if(source instanceof float[]){
			return ((float[])source).clone();
		}else{
			if(source instanceof short[]){
				return ((short[])source).clone();
			}else{
				if(source instanceof byte[]){
					return ((byte[])source).clone();
				}else{
					if(source instanceof short[])
						return ((short[])source).clone();
					else
						return ((int[])source).clone();
				}
			}
		}
	}

	public static Object[] clone_obj_array(Object[] source){
		Object[] cloned=new Object[source.length];
		for(int i=0;i<source.length;i++)
			cloned[i]=clone_obj_array(source[i]);
		return cloned;
	}

	public static Object[][] clone_obj_array(Object[][] source){
		Object[][] cloned=new Object[source.length][];
		for(int i=0;i<source.length;i++)
			cloned[i]=clone_obj_array(source[i]);
		return cloned;
	}

	public static Object[] create_stack(int length1,int length2,float value){
		Object[] temp=new Object[length1];
		for(int i=0;i<length1;i++){
			float[] temp2=new float[length2];
			if(value!=0.0f)
				Arrays.fill(temp2,value);
			temp[i]=temp2;
		}
		return temp;
	}

	public static Object[][] create_stack(int length1,int length2,int length3,float value){
		Object[][] temp=new Object[length1][];
		for(int i=0;i<length1;i++){
			temp[2]=create_stack(length2,length3,value);
		}
		return temp;
	}

	public static Object[][][] create_stack(int length1,int length2,int length3,int length4,float value){
		Object[][][] temp=new Object[length1][][];
		for(int i=0;i<length1;i++){
			temp[2]=create_stack(length2,length3,length4,value);
		}
		return temp;
	}

	public static Object transpose_image(Object source,int width,int height){
		if(source instanceof float[]){
			float[] temp=new float[width*height];
			for(int i=0;i<width;i++){
				float[] temp2=(float[])get_image_col(source,width,height,i);
				System.arraycopy(temp2,0,temp,i*height,height);
			}
			return temp;
		}else{
			if(source instanceof short[]){
				short[] temp=new short[width*height];
				for(int i=0;i<width;i++){
					short[] temp2=(short[])get_image_col(source,width,height,i);
					System.arraycopy(temp2,0,temp,i*height,height);
				}
				return temp;
			}else{
				if(source instanceof byte[]){
					byte[] temp=new byte[width*height];
					for(int i=0;i<width;i++){
						byte[] temp2=(byte[])get_image_col(source,width,height,i);
						System.arraycopy(temp2,0,temp,i*height,height);
					}
					return temp;
				}else{
					int[] temp=new int[width*height];
					for(int i=0;i<width;i++){
						int[] temp2=(int[])get_image_col(source,width,height,i);
						System.arraycopy(temp2,0,temp,i*height,height);
					}
					return temp;
				}
			}
		}
	}

	public static Object get_image_row(Object source,int width,int height,int row){
		return get_subarray(source,row*width,width);
	}

	public static Object get_image_col(Object source,int width,int height,int col){
		if(source instanceof float[]){
			float[] temp=new float[height];
			int counter=col;
			for(int i=0;i<height;i++){
				temp[i]=((float[])source)[counter];
				counter+=width;
			}
			return temp;
		}else{
			if(source instanceof short[]){
				short[] temp=new short[height];
				int counter=col;
				for(int i=0;i<height;i++){
					temp[i]=((short[])source)[counter];
					counter+=width;
				}
				return temp;
			}else{
				if(source instanceof byte[]){
					byte[] temp=new byte[height];
					int counter=col;
					for(int i=0;i<height;i++){
						temp[i]=((byte[])source)[counter];
						counter+=width;
					}
					return temp;
				}else{
					int[] temp=new int[height];
					int counter=col;
					for(int i=0;i<height;i++){
						temp[i]=((int[])source)[counter];
						counter+=width;
					}
					return temp;
				}
			}
		}
	}

	public static void set_image_col(Object source,Object dest,int width,int height,int col){
		if(dest instanceof float[]){
			float[] temp=convert_arr_float2(source);
			int counter=col;
			for(int i=0;i<height;i++){
				((float[])dest)[counter]=temp[i];
				counter+=width;
			}
		}else{
			if(dest instanceof short[]){
				short[] temp=convert_arr_short2(source);
				int counter=col;
				for(int i=0;i<height;i++){
					((short[])dest)[counter]=temp[i];
					counter+=width;
				}
			}else{
				if(dest instanceof byte[]){
					byte[] temp=convert_arr_byte2(source);
					int counter=col;
					for(int i=0;i<height;i++){
						((byte[])dest)[counter]=temp[i];
						counter+=width;
					}
				}else{
					int[] temp=convert_arr_int2(source);
					int counter=col;
					for(int i=0;i<height;i++){
						((int[])dest)[counter]=temp[i];
						counter+=width;
					}
				}
			}
		}
	}

	public static Object get_stack_col(Object[] source,int width,int height,int x,int y,int slices){
		return get_stack_col(source,width,height,x+y*width,slices);
	}

	public static Object get_stack_col(Object[] source,int width,int height,int index,int slices){
		if(source[0] instanceof float[]){
			float[] temp=new float[slices];
			for(int i=0;i<slices;i++)
				temp[i]=((float[])source[i])[index];
			return temp;
		}else{
			if(source[0] instanceof short[]){
				short[] temp=new short[slices];
				for(int i=0;i<slices;i++)
					temp[i]=((short[])source[i])[index];
				return temp;
			}else{
				if(source[0] instanceof byte[]){
					byte[] temp=new byte[slices];
					for(int i=0;i<slices;i++)
						temp[i]=((byte[])source[i])[index];
					return temp;
				}else{
					int[] temp=new int[slices];
					for(int i=0;i<slices;i++)
						temp[i]=((int[])source[i])[index];
					return temp;
				}
			}
		}
	}

	public static void set_stack_col(Object[] source,int width,int height,int x,int y,int slices,Object col){
		set_stack_col(source,width,height,x+y*width,slices,col);
	}

	public static void set_stack_col(Object[] source,int width,int height,int index,int slices,Object col){
		if(source[0] instanceof float[]){
			for(int i=0;i<slices;i++)
				((float[])source[i])[index]=((float[])col)[i];
		}else{
			if(source[0] instanceof short[]){
				for(int i=0;i<slices;i++)
					((short[])source[i])[index]=((short[])col)[i];
			}else{
				if(source[0] instanceof byte[]){
					for(int i=0;i<slices;i++)
						((byte[])source[i])[index]=((byte[])col)[i];
				}else{
					for(int i=0;i<slices;i++)
						((int[])source[i])[index]=((int[])col)[i];
				}
			}
		}
	}

	public static float get_stack_col_stat(String stat,Object[] source,int width,int height,int x,int y,int slices,float[] extras){
		Object temp=get_stack_col(source,width,height,x,y,slices);
		return jstatistics.getstatistic(stat,temp,extras);
	}

	public static float get_stack_col_stat(String stat,Object[] source,int width,int height,int index,int slices,float[] extras){
		Object temp=get_stack_col(source,width,height,index,slices);
		return jstatistics.getstatistic(stat,temp,extras);
	}

	public static float[] get_stack_proj_stat(String stat,Object[] source,int width,int height,int slices,float[] extras){
		float[] temp=new float[width*height];
		for(int i=0;i<width*height;i++)
			temp[i]=get_stack_col_stat(stat,source,width,height,i,slices,extras);
		return temp;
	}

	public static Object[] bin_stack(Object[] source,int width,int height,int binby){
		int slices=source.length;
		int newslices=(int)((float)slices/(float)binby);
		if(newslices<1){
			binby=slices;
			newslices=1;
		}
		Object[] output=new Object[newslices];
		for(int i=0;i<newslices;i++){
			Object[] tempsource=new Object[binby];
			for(int j=0;j<binby;j++){
				tempsource[j]=source[j+i*binby];
			}
			output[i]=get_stack_proj_stat("Avg",tempsource,width,height,binby,null);
		}
		return output;
	}

	public static Object[] bin_stack_stat(String stat,Object[] source,int width,int height,int start,int end,int binby,float[] extras){
		int slices=end-start+1;
		int newslices=(int)((float)slices/(float)binby);
		Object[] output=new Object[newslices];
		for(int i=0;i<newslices;i++){
			Object[] tempsource=new Object[binby];
			for(int j=0;j<binby;j++){
				tempsource[j]=source[j+i*binby+start];
			}
			output[i]=get_stack_proj_stat(stat,tempsource,width,height,binby,extras);
		}
		return output;
	}

	public static float[] pad_2D(float[] image,int width,int height,int newwidth,int newheight,int padindex){
		// here we pad an image to newwidth,newheight with 0, the avg, or the
		// edge avg
		// if newwidth or newheight is less than width or height, we crop
		float[] newimg=new float[newwidth*newheight];
		float avg=0.0f;
		if(padindex==1)
			avg=jstatistics.getstatistic("Avg",image,null);
		if(padindex==2){
			float[] edge=new float[height];
			for(int i=1;i<=height;i++)
				edge[i-1]=image[i*width-1];
			avg=jstatistics.getstatistic("Avg",edge,null);
		}
		int minheight=Math.min(height,newheight);
		int minwidth=Math.min(width,newwidth);
		for(int i=0;i<minheight;i++){
			int offset=i*newwidth;
			System.arraycopy(image,i*width,newimg,offset,minwidth);
			for(int j=width;j<newwidth;j++)
				newimg[j+offset]=avg;
		}
		float newavg=avg;
		if(padindex==2){
			float[] edge=new float[width];
			System.arraycopy(image,(height-1)*width,edge,0,width);
			newavg=jstatistics.getstatistic("Avg",edge,null);
		}
		for(int i=height;i<newheight;i++){
			for(int j=0;j<width;j++){
				newimg[j+i*newwidth]=newavg;
			}
		}
		float corneravg=0.5f*(avg+newavg);
		for(int i=height;i<newheight;i++){
			for(int j=width;j<newwidth;j++){
				newimg[j+i*newwidth]=corneravg;
			}
		}
		return newimg;
	}
	
	public static float get3DSliceStat(Object[] is,int frame,int slice,int channel,int frames,int slices,int channels,String stat){
		Object pixels=get3DSlice(is,frame,slice,channel,frames,slices,channels);
		return jstatistics.getstatistic(stat,pixels,null);
	}

	public static Object get3DSlice(Object[] is,int frame,int slice,int channel,int frames,int slices,int channels){
		return is[channel+slice*channels+frame*channels*slices];
	}
	
	public static Object get3DSliceInterp(Object[] is,int frame,float slice,int channel,int frames,int slices,int channels){
		int prev=(int)slice;
		Object sl1=get3DSlice(is,frame,prev,channel,frames,slices,channels);
		if(prev==slice) return sl1;
		int next=prev+1;
		Object sl2=get3DSlice(is,frame,next,channel,frames,slices,channels);
		return interpolation.interpz(sl1,sl2,get_array_length(sl1),slice-prev);
	}

	public static void set3DSlice(Object[] is,Object pixels,int frame,int slice,int channel,int frames,int slices,int channels){
		is[channel+slice*channels+frame*channels*slices]=pixels;
	}

	public static Object[] get3DTSeries(Object[] is,int slice,int channel,int frames,int slices,int channels){
		Object[] temp=new Object[frames];
		for(int i=0;i<frames;i++){
			temp[i]=get3DSlice(is,i,slice,channel,frames,slices,channels);
		}
		return temp;
	}

	public static void set3DTSeries(Object[] source,Object[] is,int slice,int channel,int frames,int slices,int channels){
		for(int i=0;i<frames;i++){
			set3DSlice(is,source[i],i,slice,channel,frames,slices,channels);
		}
	}

	public static Object[] get3DZSeries(Object[] is,int frame,int channel,int frames,int slices,int channels){
		Object[] temp=new Object[slices];
		for(int i=0;i<slices;i++){
			temp[i]=get3DSlice(is,frame,i,channel,frames,slices,channels);
		}
		return temp;
	}

	public static void set3DZSeries(Object[] source,Object[] is,int frame,int channel,int frames,int slices,int channels){
		for(int i=0;i<slices;i++){
			set3DSlice(is,source[i],frame,i,channel,frames,slices,channels);
		}
	}

	public static Object[] get3DCSeries(Object[] is,int frame,int slice,int frames,int slices,int channels){
		Object[] temp=new Object[channels];
		for(int i=0;i<channels;i++){
			temp[i]=get3DSlice(is,frame,slice,i,frames,slices,channels);
		}
		return temp;
	}

	public static void set3DCSeries(Object[] source,Object[] is,int frame,int slice,int frames,int slices,int channels){
		for(int i=0;i<channels;i++){
			set3DSlice(is,source[i],frame,slice,i,frames,slices,channels);
		}
	}

	public static float[] get3DProjZStat(Object[] is,int frame,int channel,int frames,int slices,int channels,String stat,float[] extras){
		Object[] zseries=get3DZSeries(is,frame,channel,frames,slices,channels);
		int len=get_array_length(zseries[0]);
		return get_stack_proj_stat(stat,zseries,len,1,slices,extras);
	}

	public static Object get3DProjZStat2(Object[] is,int frame,int channel,int frames,int slices,int channels,String stat,float[] extras){
		// here we convert back to original data type
		int typeindex=get_array_type(is[0]);
		float[] proj=get3DProjZStat(is,frame,channel,frames,slices,channels,stat,extras);
		return convert_array(proj,typeindex);
	}

	public static float[] get3DProjTStat(Object[] is,int slice,int channel,int frames,int slices,int channels,String stat,float[] extras){
		Object[] tseries=get3DTSeries(is,slice,channel,frames,slices,channels);
		int len=get_array_length(tseries[0]);
		return get_stack_proj_stat(stat,tseries,len,1,frames,extras);
	}

}
