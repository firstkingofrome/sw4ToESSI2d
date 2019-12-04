#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:48:57 2019
@author: eeckert
"""

import h5py
import re
import pandas
import argparse
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import interpolate #used for its linear interpolation function

def dframeToDict(dFrame):
    #takes the input dataframe, unfucks it, and returns a dictionary
    dFrame = list(dFrame.iterrows())
    return {i[1].to_list()[0] : i[1].to_list()[1] for i in dFrame}

def roundToBase(number,base):
    return np.int32(base*np.round(number/base))

"""
This code is probably fine, the problems must be a result of a mismatch between the mesh resolution and the data

"""
class interpolateSW4():
    def __init__(self,sw4essiout_filename,parameterFile,az):
        #load the parameter file
        self.pf = dframeToDict(pd.read_csv(parameterFile))
        #check that the azimuth is valid
        assert (int(az)==90 or int(az)==0)
        if(az==90):
            #I am using the sw4 x plane direction |
            self.directivity = 'sw4_x'
        else:
            #I am using the sw4 y plane direction --
            self.directivity = 'sw4_y'
        #verify that a valid interpolation method has been selected
        assert (self.pf["interpolate"]==None or self.pf["interpolate"]=="linear" or self.pf["interpolate"]=="cubic" or self.pf["interpolate"]=="quintic")
        #now deal with the input parameter file informations
        self.sw4essiout = h5py.File(sw4essiout_filename, 'r')
        self.h  = self.sw4essiout['ESSI xyz grid spacing'][0]
        self.x0 = self.sw4essiout['ESSI xyz origin'][0]
        self.y0 = self.sw4essiout['ESSI xyz origin'][1]
        self.z0 = self.sw4essiout['ESSI xyz origin'][2]
        self.t0 = self.sw4essiout['time start'][0]
        self.npts = self.sw4essiout['vel_0 ijk layout'].shape[0]
        self.dt = self.sw4essiout['timestep'][0]
        self.t1 = self.dt*(self.npts-1)
        self.nt = self.sw4essiout['vel_0 ijk layout'].shape[0]
        self.time = np.linspace(self.t0, self.t1, self.npts)
        if(self.directivity == 'sw4_x'):
            #load the input data
            self.plane = int(int(self.pf['y0'])/self.h)
        else:
            #load the input data
            self.plane = int(int(self.pf['x0'])/self.h)
        #now save the plane for this orientation
        if(self.directivity == 'sw4_x'):
            #load the input data
            #n_horrizontal = sw4essiout['vel_0 ijk layout'].shape[1]
            #take all timesteps, and the plane as specified in the input file
            #data_horrizontal etc. is all in ESSI COORDINATES!
            self.data_horrizontal = self.sw4essiout['vel_0 ijk layout'][:,:,self.plane,:]
            #set out of plane data to 0
            self.data_out_of_plane = np.zeros(shape=self.sw4essiout['vel_1 ijk layout'][:,:,self.plane,:].shape)
            self.data_vertical = -1.0*self.sw4essiout['vel_2 ijk layout'][:,:,self.plane,:]
        else:
            #n_horrizontal = sw4essiout['vel_1 ijk layout'].shape[0]
            self.data_horrizontal = self.sw4essiout['vel_1 ijk layout'][:,self.plane,:,:]
            #set out of plane data to 0
            self.data_out_of_plane = np.zeros(shape=self.sw4essiout['vel_0 ijk layout'][:,self.plane,:,:].shape)
            self.data_vertical = -1.0*self.sw4essiout['vel_2 ijk layout'][:,self.plane,:,:]
    def interpolateToTargetResolution(self):
        #save the horrizontal and veritical interpolation formulas in a list where index corresponds to time step
        self.hInterpolator = []
        self.vInterpolator = []
        #generate the coordinate basis so that the interpolator knows where this information is located spacially
        hMax,vMax = (self.data_horrizontal.shape[1])*self.h,(self.data_horrizontal.shape[2])*self.h
        xGrid,zGrid = np.arange(0,hMax,self.h),np.arange(0,vMax,self.h)
        #create interpolateor for the horrizontal data using this method and deal with every time step
        for t in (range(self.data_horrizontal.shape[0])):
            self.hInterpolator.append(interpolate.interp2d(xGrid, zGrid, self.data_horrizontal[t,:,:], kind=self.pf["interpolate"]))
        #now do the same for the veritcal data--The geometery of the problem is the same so I just have to change to data vertical
        for t in (range(self.data_vertical.shape[0])):
            self.vInterpolator.append(interpolate.interp2d(xGrid, zGrid, self.data_vertical[t,:,:], kind=self.pf["interpolate"]))

    def evalAtPoint(self,sw4x,sw4z):
        #creates a time series for the specified point using the specified interpolation function
        #this is then returned to the caller
        shape = (self.data_horrizontal.shape[0])
        horrizontalPlane,verticalPlane = np.zeros(shape=shape,dtype=np.float32),np.zeros(shape=shape,dtype=np.float32)
        #loop over everytime step and save the interpolated data for that time step
        for t in range(shape):
            horrizontalPlane[t] = self.hInterpolator[t](sw4x,sw4z)
            verticalPlane[t] = self.vInterpolator[t](sw4x,sw4z)
        return horrizontalPlane,verticalPlane

print("Modified!")

parameterFile="get2dSlice.csv"
sw4essiout_filename="Location17faultParallel.cycle=00000.essi"
pointsRequested="Comparison_points.csv"
pointsOutput = "outputPoints"
#load the csv file  
points = pandas.read_csv(pointsRequested)
#get some important constant values
sw4Xstart,sw4Ystart,sw4Zstart = float(points["sw4Xstart"][0]),float(points["sw4Ystart"][0]),float(points["sw4Zstart"][0])
essixStart,essiyStart,essizStart = float(points["essiXstart"][0]),float(points["essiYstart"][0]),float(points["essiZstart"][0])
#load the sw4 output an prepare the interpolators
sw4Data = interpolateSW4(sw4essiout_filename,parameterFile,int(points["azimuth"][0]))
sw4Data.interpolateToTargetResolution()
#save some constant values
directivity = sw4Data.directivity
for index, row in points.iterrows():
    print(row)
    #translate into xyz coordinates
    ### THE ESSI CORDINATES ARE XYZ WHILE SW4 IS YXZ
    OSx,OSy,OSz = row["x"],row["y"],row["z"]
    #round all of this to the nearest point as defined by the spacing
    print("preparing refrence point " + row["point"])
    print("saving x y and z velocities for  exssi point " + "%i,%i,%i" % (OSx,OSy,OSz))
    """
    THE LINE BELOW THIS IS WHERE I CONVERT FROM ESSI --> SW4
    """
    if(directivity == 'sw4_x'):
        sw4x,sw4y,sw4z = (sw4Xstart+(OSy-essiyStart)),(sw4Ystart+(OSx-essixStart)),(sw4Zstart+(essizStart-OSz))
        print("This corresponds to sw4 point %f,%f"  %(sw4x,sw4z))
        print('\n')
        #save xyz as csv files for that point
        #save ESSI x y and z!
        fname = "_sw4_"+row["point"]+".csv"
        #sample the interpolation functions to get the data
        horrizontal,vertical = sw4Data.evalAtPoint(sw4x,sw4z)
        np.savetxt(fname='x'+"_vel"+fname,X=horrizontal)
        #save y
        #y = sw4essiout['vel_0 ijk layout'][:,sw4x,sw4y,sw4z]
        #np.savetxt(fname='y'+fname,X=y)
        #save z
        np.savetxt(fname='z'+"_vel"+fname,X=vertical)

    else:
        sw4x,sw4y,sw4z = (sw4Xstart+(OSy-essiyStart)),(sw4Ystart+(OSx-essixStart)),(sw4Zstart+(essizStart-OSz))
        print("This corresponds to sw4 point %f,%f,%f"  %(sw4x,sw4y,sw4z))
        print('\n')
        #save xyz as csv files for that point
        #save ESSI x y and z!
        fname = "_sw4_"+row["point"]+".csv"
        horrizontal,vertical = sw4Data.evalAtPoint(sw4y,sw4z)
        np.savetxt(fname='x'+"_vel"+fname,X=horrizontal)
        #save y
        #y = sw4essiout['vel_0 ijk layout'][:,sw4x,sw4y,sw4z]
        #np.savetxt(fname='y'+fname,X=y)
        #save z
        np.savetxt(fname='z'+"_vel"+fname,X=vertical)



"""
if __name__ == '__main__':
    #initialize arguments
    print("starting get reference points!")
    parser = argparse.ArgumentParser()
    parser.add_argument('--sw4HDF5Out',help="current hdf5 output for sw4", nargs="+", type=str)
    parser.add_argument('--refPoints',help="csv file for input reference points", nargs="+", type=str)
    args = parser.parse_args()   
    main(args)
"""


"""
quick plot to check:
import matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(time, z_acc,color='r',label="computed acceleration")
#now differentiate the whole thing
all_data = -1.0*sw4essiout['vel_2 ijk layout'][:]
all_data_acc = np.gradient(all_data[:,:,:,:],dt,axis=0)
ax.plot(time, all_data_acc[:,int(sw4x),int(sw4y),int(sw4z)],color='b',label="computed acceleration for all points (axis=0)")
ax.legend(loc='upper left')
fig.savefig("testData.png")
"""
