/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.algutils;
import jalgs.gui_interface;
import jalgs.jstatistics;

import java.awt.Polygon;

public class findblobs3D{
	// this class finds contiguous blobs using a flood fill mechanism
	// these blobs can then be sorted using area and circularity
	// object images are typically stored as floating bit images with contiguous
	// blobs having the same index
	public int width,height,depth,nobjects;
	public findblobs3 fb;
	int maxStackSize=500; // will be increased as needed
	int[][] stack=new int[maxStackSize][3];
	int stackSize;
	public gui_interface gui;

	public findblobs3D(int width1,int height1,int depth1){
		width=width1;
		height=height1;
		depth=depth1;
		fb=new findblobs3(width,height);
		nobjects=0;
		gui=null;
	}
	
	public findblobs3D(int width1,int height1,int depth1,gui_interface gui){
		width=width1;
		height=height1;
		depth=depth1;
		fb=new findblobs3(width,height);
		nobjects=0;
		this.gui=gui;
	}
	
	public void set_objects(float[][] objects){
		nobjects=(int)maxarray(objects);
	}
	
	public float[][] dofindblobs(Object[] data1){
		if(data1[0] instanceof byte[]){
			byte[][] temp=new byte[data1.length][];
			for(int i=0;i<data1.length;i++) temp[i]=(byte[])data1[i];
			return dofindblobs(temp);
		} else {
			float[][] temp=new float[data1.length][];
			for(int i=0;i<temp.length;i++) temp[i]=(float[])data1[i];
			return dofindblobs(temp,0.5f);
		}
	}

	public float[][] dofindblobs(byte[][] data1){
		float[][] temp=new float[depth][];
		int[] sliceblobs=new int[depth];
		for(int i=0;i<depth;i++){
			temp[i]=fb.dofindblobs(data1[i]);
			sliceblobs[i]=fb.nobjects;
		}
		float[][] objects=new float[depth][width*height];
		// now go through and assemble the 2D objects into 3D objects
		int id=0;
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(temp[i][k+j*width]>0.0f){
						id++;
						fill3D(temp,objects,sliceblobs,id,k,j,i);
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		nobjects=id;
		return objects;
	}

	public float[][] dofindblobs(float[][] data1,float thresh){
		float[][] temp=new float[data1.length][];
		int[] sliceblobs=new int[depth];
		for(int i=0;i<data1.length;i++){
			temp[i]=fb.dofindblobs(data1[i],thresh);
			sliceblobs[i]=fb.nobjects;
		}
		float[][] objects=new float[depth][width*height];
		// now go through and assemble the 2D objects into 3D objects
		int id=0;
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(temp[i][k+j*width]>0.0f){
						id++;
						fill3D(temp,objects,sliceblobs,id,k,j,i);
					}
				}
			}
		}
		nobjects=id;
		return objects;
	}
	
	public byte[][] tobinary(float[][] objects){
		return tobinary(objects,false);
	}

	public byte[][] tobinary(float[][] objects,boolean separate){
		byte[][] out=new byte[depth][width*height];
		if(separate) separateobjects(objects);
		for(int i=0;i<depth;i++){
			for(int j=0;j<width*height;j++)
				if(objects[i][j]>0.0f)
					out[i][j]=(byte)255;
		}
		return out;
	}

	public void clear_edges(float[][] objects,boolean renumber){
		int totlength=width*height+width*depth+height*depth;
		int counter=0;
		for(int i=0;i<width*height;i++){
			// get the top and bottom
			if(objects[0][i]>0.0f){
				delete_object(objects,objects[0][i]);
			}
			if(objects[depth-1][i]>0.0f){
				delete_object(objects,objects[depth-1][i]);
			}
			counter++;
			if(gui!=null) gui.showProgress(counter,totlength);
		}
		for(int k=1;k<(depth-1);k++){
			// now the sides
			for(int i=0;i<width;i++){
				if(objects[k][i]>0.0f){
					delete_object(objects,objects[k][i]);
				}
				if(objects[k][i+(height-1)*width]>0.0f){
					delete_object(objects,objects[k][i+(height-1)*width]);
				}
				counter++;
				if(gui!=null) gui.showProgress(counter,totlength);
			}
			for(int i=0;i<height;i++){
				if(objects[k][i*width]>0.0f){
					delete_object(objects,objects[k][i*width]);
				}
				if(objects[k][i*width+width-1]>0.0f){
					delete_object(objects,objects[k][i*width+width-1]);
				}
				counter++;
				if(gui!=null) gui.showProgress(counter,totlength);
			}
		}
		if(renumber) renumber_objects(objects);
	}
	
	public void separateobjects(float[][] objects){
		//here non-zero pixels neighboring another object are turned black
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					float val=cloned[i][k+j*width];
					if(val>0.0f){
						float[] neighbors=getNeighbors(cloned,k,j,i);
						for(int l=0;l<neighbors.length;l++){
							if(neighbors[l]>0.0f && neighbors[l]!=val){
								objects[i][j*width+k]=0.0f;
								break;
							}
						}
					}
				}
			}
		}
		/*byte[][] temp=tobinary(objects,false);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);*/
	}
	
	public void erodeobjects(float[][] objects,int threshold){
		//pixels neighbored in 3D by space or another object are deleted
		//objects are allowed to separate (become multiple new objects)
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					float val=cloned[i][k+j*width];
					float[] neighbors=getNeighbors(cloned,k,j,i);
					int count=0;
					for(int l=0;l<neighbors.length;l++) if(neighbors[l]!=val) count++;
					if(count>threshold) objects[i][j*width+k]=0.0f;
				}
			}
		}
		byte[][] temp=tobinary(objects,false);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);
	}
	
	public void erodeobjects2(float[][] objects,int threshold){
		//pixels neighbored in 3D by space or another object are deleted
		//objects are not allowed to separate
		//they are not renumbered either
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					float val=cloned[i][k+j*width];
					float[] neighbors=getNeighbors(cloned,k,j,i);
					int count=0;
					for(int l=0;l<neighbors.length;l++) if(neighbors[l]!=val) count++;
					if(count>threshold) objects[i][j*width+k]=0.0f;
				}
			}
		}
		//byte[][] temp=tobinary(objects,false);
		//float[][] temp2=dofindblobs(temp);
		//copy_objects(temp2,objects);
	}
	
	public void dilateobjects(float[][] objects){
		//pixels neighbored in 3D by one object are added to that object
		//objects are not allowed to merge;
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					if(cloned[i][k+j*width]==0.0f){
    					float[] neighbors=getNeighbors(cloned,k,j,i);
    					float val=0.0f; boolean found=false; boolean found2=false;
    					for(int l=0;l<neighbors.length;l++){
    						if(neighbors[l]>0.0f){
    							if(!found){
    								found=true;
    								val=neighbors[l];
    							} else {
    								if(neighbors[l]!=val){
    									found2=true; break;
    								}
    							}
    						}
    					}
    					if(found && !found2) objects[i][j*width+k]=val;
					}
				}
			}
		}
		byte[][] temp=tobinary(objects);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);
	}
	
	public void dilateobjects2(float[][] objects){
		//pixels neighbored in 3D by one object are added to that object
		//here objects are allowed to merge
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					if(cloned[i][k+j*width]==0.0f){
    					float[] neighbors=getNeighbors(cloned,k,j,i);
    					for(int l=0;l<neighbors.length;l++){
    						if(neighbors[l]>0.0f){
    							objects[i][j*width+k]=1.0f;
    							break;
    						}
    					}
					}
				}
			}
		}
		byte[][] temp=tobinary(objects,false);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);
	}
	
	private void copy_objects(float[][] source,float[][] dest){
		for(int i=0;i<source.length;i++){
			System.arraycopy(source[i],0,dest[i],0,source[i].length);
		}
	}

	public void delete_object(float[][] objects,float id){
		nobjects--;
		for(int k=0;k<depth;k++){
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					if(objects[k][j+i*width]==id){
						objects[k][j+i*width]=0.0f;
					}
				}
			}
		}
	}

	private float maxarray(float[] input){
		float max=input[0];
		for(int i=1;i<input.length;i++){
			if(input[i]>max){
				max=input[i];
			}
		}
		return max;
	}

	private float maxarray(float[][] input){
		float max=maxarray(input[0]);
		for(int i=1;i<input.length;i++){
			float tm=maxarray(input[i]);
			if(tm>max)
				max=tm;
		}
		return max;
	}

	public int get_nblobs(float[][] objects){
		return (int)maxarray(objects);
	}

	public void renumber_objects(float[][] objects){
		int maxid=(int)maxarray(objects);
		boolean[] occupancy=new boolean[maxid+1];
		for(int j=0;j<depth;j++){
			for(int i=0;i<width*height;i++){
				if(objects[j][i]>0.0f){
					occupancy[(int)objects[j][i]]=true;
				}
			}
		}
		int counter=0;
		float[] newids=new float[maxid+1];
		for(int i=1;i<=maxid;i++){
			if(occupancy[i]){
				counter++;
				newids[i]=counter;
			}
		}
		nobjects=counter;
		for(int j=0;j<depth;j++){
			for(int i=0;i<width*height;i++){
				if(objects[j][i]>0.0f){
					objects[j][i]=newids[(int)objects[j][i]];
				}
			}
		}
	}

	public float[] getNeighbors(float[][] objects,int x,int y,int z){
		//this excludes the current pixel
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		if(z==0||z>=(depth-1)){
			return null;
		}
		float[] temp=new float[26];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[z-1][temp2];
		temp2++;
		temp[1]=objects[z-1][temp2];
		temp2++;
		temp[2]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[3]=objects[z-1][temp2];
		temp2++;
		temp[4]=objects[z-1][temp2];
		temp2++;
		temp[5]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[6]=objects[z-1][temp2];
		temp2++;
		temp[7]=objects[z-1][temp2];
		temp2++;
		temp[8]=objects[z-1][temp2];
		temp2=x-1+(y-1)*width;
		temp[9]=objects[z][temp2];
		temp2++;
		temp[10]=objects[z][temp2];
		temp2++;
		temp[11]=objects[z][temp2];
		temp2+=(width-2);
		temp[12]=objects[z][temp2];
		temp2+=2;
		temp[13]=objects[z][temp2];
		temp2+=(width-2);
		temp[14]=objects[z][temp2];
		temp2++;
		temp[15]=objects[z][temp2];
		temp2++;
		temp[16]=objects[z][temp2];
		temp2=x-1+(y-1)*width;
		temp[17]=objects[z+1][temp2];
		temp2++;
		temp[18]=objects[z+1][temp2];
		temp2++;
		temp[19]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[20]=objects[z+1][temp2];
		temp2++;
		temp[21]=objects[z+1][temp2];
		temp2++;
		temp[22]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[23]=objects[z+1][temp2];
		temp2++;
		temp[24]=objects[z+1][temp2];
		temp2++;
		temp[25]=objects[z+1][temp2];
		return temp;
	}

	public float[] getNeighbors2(float[][] objects,int x,int y,int z){
		//here we include the center pixel
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		if(z==0||z>=(depth-1)){
			return null;
		}
		float[] temp=new float[27];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[z-1][temp2];
		temp2++;
		temp[1]=objects[z-1][temp2];
		temp2++;
		temp[2]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[3]=objects[z-1][temp2];
		temp2++;
		temp[4]=objects[z-1][temp2];
		temp2++;
		temp[5]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[6]=objects[z-1][temp2];
		temp2++;
		temp[7]=objects[z-1][temp2];
		temp2++;
		temp[8]=objects[z-1][temp2];
		temp2=x-1+(y-1)*width;
		temp[9]=objects[z][temp2];
		temp2++;
		temp[10]=objects[z][temp2];
		temp2++;
		temp[11]=objects[z][temp2];
		temp2+=(width-2);
		temp[12]=objects[z][temp2];
		temp2++;
		temp[13]=objects[z][temp2];
		temp2++;
		temp[14]=objects[z][temp2];
		temp2+=(width-2);
		temp[15]=objects[z][temp2];
		temp2++;
		temp[16]=objects[z][temp2];
		temp2++;
		temp[17]=objects[z][temp2];
		temp2=x-1+(y-1)*width;
		temp[18]=objects[z+1][temp2];
		temp2++;
		temp[19]=objects[z+1][temp2];
		temp2++;
		temp[20]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[21]=objects[z+1][temp2];
		temp2++;
		temp[22]=objects[z+1][temp2];
		temp2++;
		temp[23]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[24]=objects[z+1][temp2];
		temp2++;
		temp[25]=objects[z+1][temp2];
		temp2++;
		temp[26]=objects[z+1][temp2];
		return temp;
	}

	// Does a 4-connected 3D flood fill
	// the data input stack has been flood filled already with the 2D algorithm
	// the object gets deleted from data
	public boolean fill3D(float[][] data,float[][] counter,int[] sliceblobs,int id,int x,int y,int z){
		stackSize=0;
		push(x,y,z);
		while(true){
			int[] cds=pop();
			x=cds[0];
			if(x<0)
				return true;
			y=cds[1];
			z=cds[2];
			if(data[z][x+y*width]==0.0f)
				continue;
			int[] lims2d=fb.getfilllimits(data[z],(int)data[z][x+y*width],x,y);
			//if(lims2d[0]<0) lims2d[0]=0;  if(lims2d[1]>=width) lims2d[1]=(width-1);
			//if(lims2d[2]<0) lims2d[2]=0;  if(lims2d[3]>=height) lims2d[3]=(height-1);
			fillRegion(data,counter,id,lims2d,z,(int)data[z][x+y*width]); // fill scan-region
			if(z>0){
				boolean[] prevselected=new boolean[sliceblobs[z-1]];
				for(int i=lims2d[2];i<=lims2d[3];i++){ // find unique scan-regions above this one
					for(int j=lims2d[0];j<=lims2d[1];j++){
						if(counter[z][j+i*width]==id){
							int tempid=(int)data[z-1][j+i*width];
							if(tempid>0 && !prevselected[tempid-1]){
								prevselected[tempid-1]=true;
								push(j,i,z-1);
							}
						}
					}
				}
			}
			if(z<(depth-1)){
				boolean[] prevselected=new boolean[sliceblobs[z+1]];
				for(int i=lims2d[2];i<=lims2d[3];i++){ // find unique scan-regions below this one
					for(int j=lims2d[0];j<=lims2d[1];j++){
						if(counter[z][j+i*width]==id){
							int tempid=(int)data[z+1][j+i*width];
							if(tempid>0&&!prevselected[tempid-1]){
								prevselected[tempid-1]=true;
								push(j,i,z+1);
							}
						}
					}
				}
			}
		}
	}

	final void fillRegion(float[][] data,float[][] counter,int id,int[] lims2d,int z,int oldid){
		for(int i=lims2d[2];i<=lims2d[3];i++){
			for(int j=lims2d[0];j<=lims2d[1];j++){
				if(data[z][j+i*width]==oldid){
					data[z][j+i*width]=0.0f;
					counter[z][j+i*width]=id;
				}
			}
		}
	}

	final void push(int x,int y,int z){
		stackSize++;
		if(stackSize==maxStackSize){
			int[][] newStack=new int[maxStackSize*2][];
			System.arraycopy(stack,0,newStack,0,maxStackSize);
			stack=newStack;
			maxStackSize*=2;
		}
		stack[stackSize-1]=new int[]{x,y,z};
	}

	final int[] pop(){
		if(stackSize<1)
			return new int[]{-1,-1,-1};
		else{
			int[] vals=stack[stackSize-1];
			stackSize--;
			return vals;
		}
	}
	
	public int[][] getallfilllimits(float[][] objects){
		//int nobjects=(int)maxarray(objects);
		//the limits are lowerx,upperx,lowery,uppery,lowerz,upperz
		int[][] lims=new int[nobjects][6];
		for(int i=0;i<nobjects;i++){
			lims[i][0]=width-1; lims[i][1]=0;
			lims[i][2]=height-1; lims[i][3]=0;
			lims[i][4]=depth-1; lims[i][5]=0;
		}
		for(int k=0;k<depth;k++){
    		for(int i=0;i<height;i++){
    			for(int j=0;j<width;j++){
    				if(objects[k][j+i*width]>0.0f){
    					int id=(int)objects[k][j+i*width]-1;
    					if(j<lims[id][0]) lims[id][0]=j;
    					if(j>lims[id][1]) lims[id][1]=j;
    					if(i<lims[id][2]) lims[id][2]=i;
    					if(i>lims[id][3]) lims[id][3]=i;
    					if(k<lims[id][4]) lims[id][4]=k;
    					if(k>lims[id][5]) lims[id][5]=k;
    				}
    			}
    		}
		}
		for(int i=0;i<nobjects;i++){
			if(lims[i][1]<lims[i][0]){
				//here we never found an object, so reset full limits
				lims[i][0]=0; lims[i][1]=width-1;
				lims[i][2]=0; lims[i][3]=height-1;
				lims[i][4]=depth-1; lims[i][5]=0;
			}
		}
		return lims;
	}
	
	/********************************************************
	 * gets all of of the object stats
	 * @param objects
	 * @param measurement: the intensity image to measure
	 * @param lims
	 * @param stat
	 * @return
	 */
	public float[] get_all_object_stats(float[][] objects,Object[] measurement,int[][] lims,String stat){
		float[] stats=new float[nobjects];
		int[] areas=get_areas(objects);
		for(int i=0;i<nobjects;i++){
			float[] temp=new float[areas[i]];
			int counter=0;
			if(measurement[0] instanceof float[]){
    			for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							temp[counter]=((float[])measurement[j])[xyindex];
    							counter++;
    						}
    					}
    				}
    			}
			} else if(measurement[0] instanceof short[]){
				for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							temp[counter]=((short[])measurement[j])[xyindex]&0xffff;
    							counter++;
    						}
    					}
    				}
    			}
			} else if(measurement[0] instanceof byte[]){
				for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							temp[counter]=((byte[])measurement[j])[xyindex]&0xff;
    							counter++;
    						}
    					}
    				}
    			}
			}
			stats[i]=jstatistics.getstatistic(stat,temp,null);
			if(gui!=null) gui.showProgress(i,nobjects);
		}
		return stats;
	}
	
	public int[] get_areas(float[][] objects){
		int nblobs=get_nblobs(objects);
		int[] hist=new int[nblobs];
		for(int j=0;j<depth;j++){
    		for(int i=0;i<width*height;i++){
    			if(objects[j][i]>0.0f){
    				hist[(int)objects[j][i]-1]++;
    			}
    		}
		}
		return hist;
	}
	
	public int[] get_areas(float[][] objects,float[][] lims){
		int nblobs=get_nblobs(objects);
		int[] hist=new int[nblobs];
		for(int j=0;j<depth;j++){
    		for(int i=0;i<width*height;i++){
    			if(objects[j][i]>0.0f){
    				hist[(int)objects[j][i]-1]++;
    			}
    		}
		}
		return hist;
	}
	
	public void filter_area(float[][] objects,int[] arealims,boolean renumber){
		int[] hist=get_areas(objects);
		for(int i=0;i<hist.length;i++){
			if(hist[i]<arealims[0]||hist[i]>arealims[1]){
				delete_object(objects,i+1);
			}
			if(gui!=null) gui.showProgress(i,hist.length);
		}
		if(renumber)
			renumber_objects(objects);
	}
	
	public float[][] getCentroidsAreas(float[][] objects){
		int nblobs=get_nblobs(objects);
		float[][] centroids=new float[nblobs][4];
		for(int i=0;i<objects.length;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					int index=(int)objects[i][k+j*width];
					if(index>0){
						centroids[index-1][0]+=k;
						centroids[index-1][1]+=j;
						centroids[index-1][2]+=i;
						centroids[index-1][3]+=1.0f;
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		for(int i=0;i<nblobs;i++){
			centroids[i][0]/=centroids[i][3];
			centroids[i][1]/=centroids[i][3];
			centroids[i][2]/=centroids[i][3];
		}
		return centroids;
	}
	
	public float[][] getCentroidsAreasAvgs(float[][] objects,Object[] measurement){
		int nblobs=get_nblobs(objects);
		float[][] centroids=new float[nblobs][5];
		for(int i=0;i<objects.length;i++){
			float[] tempmeas=algutils.convert_arr_float2(measurement[i]);
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					int index=(int)objects[i][k+j*width];
					if(index>0){
						centroids[index-1][0]+=k;
						centroids[index-1][1]+=j;
						centroids[index-1][2]+=i;
						centroids[index-1][3]+=1.0f;
						centroids[index-1][4]+=tempmeas[k+j*width];
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		for(int i=0;i<nblobs;i++){
			centroids[i][0]/=centroids[i][3];
			centroids[i][1]/=centroids[i][3];
			centroids[i][2]/=centroids[i][3];
			centroids[i][4]/=centroids[i][3];
		}
		return centroids;
	}
	
	public Object[] get_object_outline(float[][] objects,int id){
		//this is incorrect--need to look for multiple polygons per slice
		Polygon[] outlines=new Polygon[objects.length];
		int[] pos=new int[objects.length];
		int counter=0;
		for(int i=0;i<objects.length;i++){
			Polygon temp=fb.get_object_outline(objects[i],id);
			if(temp!=null){
				outlines[counter]=temp;
				pos[counter]=i;
				counter++;
			}
		}
		Polygon[] outlines2=new Polygon[counter];
		System.arraycopy(outlines,0,outlines2,0,counter);
		int[] pos2=new int[counter];
		System.arraycopy(pos,0,pos2,0,counter);
		return new Object[]{outlines2,pos2};
	}

}
