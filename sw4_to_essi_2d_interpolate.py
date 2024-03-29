#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  My newley cleaned up and legible version of the interpolation codes
#  
#  Copyright 2019 Eric Eckert <dukeleto@spica>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#   
import sys
import h5py
import numpy as np
import datetime
import pandas as pd
import itertools 
import scipy
from scipy import integrate
from scipy import interpolate #used for its linear interpolation function


class interpolateSW4():
    def __init__(self,sw4essiout_filename,parameterFile):
        #load the parameter file
        self.pf = parameterFile
        az = int(self.pf['azimuth'])
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
    def interpolate(self):
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
    def closeContainer(self):
        self.sw4essiout.close()

class geometeryFile():
    def __init__(self,drm_filename,parameterFile,timeSeries,nt,az):
        self.drm_file = h5py.File(drm_filename,"r+")
        #self.sw4h = sw4h
        self.pf = parameterFile
        self.x0,self.y0 = float(parameterFile["x0"]),float(parameterFile["y0"])
        self.essiXstart,self.essiYstart,self.essiZstart = float(parameterFile["essiXstart"]),float(parameterFile["essiYstart"]),float(parameterFile["essiZstart"])
        #direction
        assert (int(az)==90 or int(az)==0)
        if(az==90):
            #I am using the sw4 x plane direction |
            self.directivity = 'sw4_x'
        elif(az==0):
            #I am using the sw4 y plane direction --
            self.directivity = 'sw4_y'
        self.coordinates = self.drm_file['Coordinates'][:]
        self.n_coord = int(self.coordinates.shape[0] / 3)
        #convert the coordinates to a list
        self.coordinates = list(self.coordinates.reshape(self.coordinates.shape[0]//3,-1,3))
        #create data sets for outputs (as necessary)
        if not "Time" in self.drm_file.keys():
            self.Time_dset          = self.drm_file.create_dataset("Time", data=timeSeries)
            self.Accelerations_dset = self.drm_file.create_dataset("Accelerations", (self.n_coord*3, nt))
            self.Displacements_dset = self.drm_file.create_dataset("Displacements", (self.n_coord*3, nt))
        else:
            del self.drm_file['Time']
            self.Time_dset = self.drm_file.create_dataset("Time", data=timeSeries)
            del self.drm_file['Accelerations']
            self.Accelerations_dset = self.drm_file.create_dataset("Accelerations", (self.n_coord*3, nt))
            del self.drm_file['Displacements']
            self.Displacements_dset = self.drm_file.create_dataset("Displacements", (self.n_coord*3, nt))
        
    def generateSW4coordinates(self):
        #convert ess coordinates to there corresponding sw4 coordinates
        self.convertedCorrdinates = []
        for coordinates in self.coordinates:
            #check the orientation
            if(self.directivity=='sw4_x'):
                x,z = coordinates[0][1]-self.essiYstart+self.x0,self.essiZstart-coordinates[0][2]
                self.convertedCorrdinates.append(np.array([x,z]))
            
            elif(self.directivity=='sw4_y'):
                x,z = coordinates[0][0]-self.essiXstart+self.y0, self.essiZstart-coordinates[0][2]
                self.convertedCorrdinates.append(np.array([x,z]))

    def saveMotions(self,interpolator):
        #evaluates the motions using the sw4 interpolator
        for i in range(len(self.convertedCorrdinates)):
            #save the time series for this point saving horrizontal and vertical velocities
            v_h,v_v = interpolator.evalAtPoint(self.convertedCorrdinates[i][0],self.convertedCorrdinates[i][1])
            #save the acceleration
            a_h = np.gradient(v_h[:],interpolator.dt,axis=0)
            a_v = np.gradient(v_v [:],interpolator.dt,axis=0)
            #save the displacement
            disp_h = scipy.integrate.cumtrapz(y=v_h[:],dx=interpolator.dt,initial=0,axis=0)
            disp_v = scipy.integrate.cumtrapz(y=v_v[:],dx=interpolator.dt,initial=0,axis=0)
            #assign to ESSI data set based on the current orientation
            if(self.directivity=='sw4_x')
                self.Accelerations_dset[i*3,:] = 0.0
                self.Accelerations_dset[i*3+1,:] = a_h 
                self.Accelerations_dset[i*3+2,:] = a_v
                
                self.Displacements_dset[i*3,:] = 0.0
                self.Displacements_dset[i*3+1,:] = disp_h 
                self.Displacements_dset[i*3+2,:] = disp_v
            else:
                self.Accelerations_dset[i*3,:] = a_h 
                self.Accelerations_dset[i*3+1,:] = 0.0
                self.Accelerations_dset[i*3+2,:] = a_v
                
                self.Displacements_dset[i*3,:] = disp_h 
                self.Displacements_dset[i*3+1,:] = 0.0
                self.Displacements_dset[i*3+2,:] = disp_v
                            
            #stop
            
            break
    def closeContainer(self):
        self.drm_file.close()

"""
Check that the integration and differentiation look reasonalbe for a point
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
ax.plot(time, disp_v)
fig.savefig("disp_plot_test.png")

fig, ax = plt.subplots()
ax.plot(time, a_v)
fig.savefig("acc_plot_test.png")
"""
def dframeToDict(dFrame):
    #takes the input dataframe, unfucks it, and returns a dictionary
    dFrame = list(dFrame.iterrows())
    return {i[1].to_list()[0] : i[1].to_list()[1] for i in dFrame}

pfName = "get2dSlice.csv"
geomFileName = "DRM_input_Model24.hdf5"
sw4MotionsFileName = "Location17faultParallel.cycle=00000.essi"
#load up the parameter file
parameterFile = dframeToDict(pd.read_csv(pfName))
#load up the sw4 motions
sw4Data = interpolateSW4(sw4MotionsFileName,parameterFile)
#prepare the interpolator
sw4Data.interpolate()
#create the time series for essi
time = np.linspace(sw4Data.t0, sw4Data.t1, sw4Data.npts)
#load up the geometery file
essiGeom = geometeryFile(geomFileName,parameterFile,time,sw4Data.npts,int(parameterFile["azimuth"]))
print('# of coordinates: ', essiGeom.n_coord)
#convert all essi coordinates to sw4 coordinates
essiGeom.generateSW4coordinates()
#go through the geometery file, evaluate every point and save the resulting time series
essiGeom.saveMotions(sw4Data)
#close both containers
essiGeom.closeContainer()
sw4Data.closeContainer()


























"""

All of this stuff pertains to the old version

def dframeToDict(dFrame):
    #takes the input dataframe, unfucks it, and returns a dictionary
    dFrame = list(dFrame.iterrows())
    return {i[1].to_list()[0] : i[1].to_list()[1] for i in dFrame}
          
#stuff to insert into main

sw4essiout_filename = "Location17faultParallel.cycle=00000.essi"
drm_filename = "DRM_input_Model27.hdf5"
parameterFile = "get2dSlice.csv"
#load up the parameter file
parameterFile = dframeToDict(pd.read_csv(parameterFile))
#check that the azimuth is valid
assert (int(parameterFile['azimuth'])==90 or int(parameterFile['azimuth'])==0)
if(parameterFile['azimuth']==90):
    #I am using the sw4 x plane direction |
    directivity = 'sw4_x'
else:
    #I am using the sw4 y plane direction --
    directivity = 'sw4_y'
#verify that a valid interpolation method has been selected
assert (parameterFile["interpolate"]==None or parameterFile["interpolate"]=="linear" or parameterFile["interpolate"]=="cubic" or parameterFile["interpolate"]=="quintic")

#load up the sw4 outputs
sw4essiout = h5py.File(sw4essiout_filename, 'r')
h  = sw4essiout['ESSI xyz grid spacing'][0]
x0 = sw4essiout['ESSI xyz origin'][0]
y0 = sw4essiout['ESSI xyz origin'][1]
z0 = sw4essiout['ESSI xyz origin'][2]
t0 = sw4essiout['time start'][0]
npts = sw4essiout['vel_0 ijk layout'].shape[0]
dt = sw4essiout['timestep'][0]
t1 = dt*(npts-1)
print ('grid spacing, h: ', h)
print ('ESSI origin x0, y0, z0: ', x0, y0, z0)
print('timing, t0, dt, npts, t1: ', t0, round(dt,6), npts, round(t1,6) )
print('Shape of HDF5 data: ', sw4essiout['vel_0 ijk layout'].shape)

nt = sw4essiout['vel_0 ijk layout'].shape[0]
time = np.linspace(t0, t1, npts)

if(directivity == 'sw4_x'):
 #load the input data
 plane = int(int(parameterFile['y0'])/h)
else:
 #load the input data
 plane = int(int(parameterFile['x0'])/h)

if(directivity == 'sw4_x'):
 #load the input data
    #n_horrizontal = sw4essiout['vel_0 ijk layout'].shape[1]
    #take all timesteps, and the plane as specified in the input file
    #data_horrizontal etc. is all in ESSI COORDINATES!
    data_horrizontal = sw4essiout['vel_0 ijk layout'][:,:,plane,:]
    #set out of plane data to 0
    data_out_of_plane = np.zeros(shape=sw4essiout['vel_1 ijk layout'][:,:,plane,:].shape)
    data_vertical = -1.0*sw4essiout['vel_2 ijk layout'][:,:,plane,:]
else:

    #n_horrizontal = sw4essiout['vel_1 ijk layout'].shape[0]
    data_horrizontal = sw4essiout['vel_1 ijk layout'][:,plane,:,:]
    #set out of plane data to 0
    data_out_of_plane = np.zeros(shape=sw4essiout['vel_0 ijk layout'][:,plane,:,:].shape)
    data_vertical = -1.0*sw4essiout['vel_2 ijk layout'][:,plane,:,:]

"""
"""
#check the spacing requested and the input spacing, if they disagree interpolate the inpute data!
if(int(parameterFile["target_h"]) != int(sw4essiout["ESSI xyz grid spacing"][0])):
    #interpolate with the specified algorithm
    print("Interpolating data from " + str(sw4essiout["ESSI xyz grid spacing"][0])
     + " to " + str(parameterFile["target_h"]) )
    
    #
"""
"""
#load up the geometery file
drm_file = h5py.File(drm_filename,"a")
coordinates = drm_file['Coordinates'][:]
n_coord = int(coordinates.shape[0] / 3)
print('# of coordinates: ', n_coord)
#create data sets for outputs (as necessary)
if not "Time" in drm_file.keys():
    Time_dset          = drm_file.create_dataset("Time", data=time)
    Accelerations_dset = drm_file.create_dataset("Accelerations", (n_coord*3, nt))
    Displacements_dset = drm_file.create_dataset("Displacements", (n_coord*3, nt))
else:
    del drm_file['Time']
    Time_dset = drm_file.create_dataset("Time", data=time)
    del drm_file['Accelerations']
    Accelerations_dset = drm_file.create_dataset("Accelerations", (n_coord*3, nt))
    del drm_file['Displacements']
    Displacements_dset = drm_file.create_dataset("Displacements", (n_coord*3, nt))

#check if I need to convert the essi coordinates to a different unit
if(parameterFile["essi_units"]=="feet"):
    print("Warning the code to preform feet to meters conversion remains untested")
    #convert the input coordinates to meters--Rounding to nearest int and casting
    print("interpolating in feet requires that the interpolation option be used, now preparing to interpolate!")
    coordinates = coordinates/3.281
    precision = np.int32(parameterFile["precision"])
    #multiplying coordinates to target precision (since I really need to be working in whole integers)
    coordinates = np.int32(coordinates*precision*10)
    print("Interpolating from " + str(h) + " ft grid spacing to " + parameterFile["target_spacing"] + "ft grid spacing")
    #target h is in feet
    target_h = np.int32((np.float32(parameterFile["target_spacing"])/3.281)*np.int32(precision)*10)
    #h is all ready in meters so just  X10Xprecision
    h = h*precision*10
    #prepare to interpolate using the specified routine 
    #the *Grids are used to construct a mesh grid of all of the possible points in the space
    #data horizontal and data vertical have the same shapes
    hMax,vMax = data_horrizontal.shape[1]*h,data_horrizontal.shape[2]*h
    xGrid,zGrid = np.arange(0,hMax,h),np.arange(0,vMax,h)
    #new grid dimensions to evaluate on
    # seems to think that splines are a good idea https://link.springer.com/content/pdf/10.1007%2Fs13201-018-0807-6.pdf
    newXGrid,newZGrid= np.arange(0,hMax,target_h),np.arange(0,vMax,target_h)
    interpolatedPlaneHorrizontal = np.zeros(shape=(data_horrizontal.shape[0],newXGrid.shape[0],newZGrid.shape[0]),dtype=np.float32)
    interpolatedPlaneVertical = np.zeros(shape=(data_horrizontal.shape[0],newXGrid.shape[0],newZGrid.shape[0]),dtype=np.float32)
    #create interpolateor for the horrizontal data using this method and deal with every time step
    for t in (range(data_horrizontal.shape[0])):
        interpolator = interpolate.interp2d(xGrid, zGrid, data_horrizontal[t,:,:], kind=parameterFile["interpolate"])
        #evaluate this every where to get the new interpolated grid 
        interpolatedPlaneHorrizontal[t,:,:] = interpolator(newXGrid,newZGrid)
    #now do the same for the veritcal data--The geometery of the problem is the same so I just have to change to data vertical
    for t in (range(data_horrizontal.shape[0])):
        interpolator = interpolate.interp2d(xGrid, zGrid, data_vertical[t,:,:], kind=parameterFile["interpolate"])
        #evaluate this every where to get the new interpolated grid 
        interpolatedPlaneVertical[t,:,:] = interpolator(newXGrid,newZGrid)
    #test it by plotting what the interpolated plane looks like
    
"""
    
"""
    #test code to plot the above and verify that it functions correctly:
    import matplotlib.pyplot as plt
    #identify the time step with the maxium average value and plot there
    t=np.argmax(data_horrizontal[:,0,0]) #time step to plot at
    plt.subplot(221)
    plt.imshow(interpolatedPlaneHorrizontal[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Horriz. Interp. '+ parameterFile["interpolate"] + " when t=" + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #now do the other one
    plt.subplot(222)
    plt.imshow(data_horrizontal[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Orig. Horriz. at time=' + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #plot the vertical plane itnerpolations
    plt.subplot(223)
    plt.imshow(interpolatedPlaneVertical[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Vertical Interp. '+ parameterFile["interpolate"] + " when t=" + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #now do the other one
    plt.subplot(224)
    plt.imshow(data_vertical[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Vert. orig. at time=' + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #save it
    plt.savefig("origionalVsInterpolated.png",dpi=400)
    plt.close()
"""
"""
    #assign data_horrizontal and data_vertical to the interpolated data
    data_horrizontal,data_vertical = interpolatedPlaneHorrizontal,interpolatedPlaneVertical
    
#Otherwise see if I need to interpolate meters
#this is the same procedure as for the other data, I just dont convert to feet
elif(parameterFile["interpolate"] != None and parameterFile["essi_units"]=="meter"):
    #convert the input coordinates to meters--Rounding to nearest int and casting
    coordinates = coordinates
    precision = np.int32(10**np.int32(parameterFile["precision"]))
    parameterFile['essi_z_end'] = np.int32(np.float32(parameterFile['essi_z_end'])*precision)
    #multiplying coordinates to target precision (since I really need to be working in whole integers)
    coordinates = np.int32(coordinates*precision)
    print("Interpolating from " + str(h) + " m grid spacing to " + parameterFile["target_spacing"] + "m grid spacing")
    #target h is in feet
    target_h = np.int32((np.float32(parameterFile["target_spacing"]))*precision)
    #h is all ready in meters so just  X10Xprecision
    h = h*precision
    #prepare to interpolate using the specified routine 
    #the *Grids are used to construct a mesh grid of all of the possible points in the space
    #data horizontal and data vertical have the same shapes
    hMax,vMax = (data_horrizontal.shape[1])*h,(data_horrizontal.shape[2])*h
    xGrid,zGrid = np.arange(0,hMax,h),np.arange(0,vMax,h)
    #new grid dimensions to evaluate on
    # seems to think that splines are a good idea https://link.springer.com/content/pdf/10.1007%2Fs13201-018-0807-6.pdf
    newXGrid,newZGrid= np.arange(0,hMax,target_h),np.arange(0,vMax,target_h)
    interpolatedPlaneHorrizontal = np.zeros(shape=(data_horrizontal.shape[0],newXGrid.shape[0],newZGrid.shape[0]),dtype=np.float32)
    interpolatedPlaneVertical = np.zeros(shape=(data_horrizontal.shape[0],newXGrid.shape[0],newZGrid.shape[0]),dtype=np.float32)
    #create interpolateor for the horrizontal data using this method and deal with every time step
    for t in (range(data_horrizontal.shape[0])):
        interpolator = interpolate.interp2d(xGrid, zGrid, data_horrizontal[t,:,:], kind=parameterFile["interpolate"])
        #evaluate this every where to get the new interpolated grid 
        interpolatedPlaneHorrizontal[t,:,:] = interpolator(newXGrid,newZGrid)
    #now do the same for the veritcal data--The geometery of the problem is the same so I just have to change to data vertical
    for t in (range(data_horrizontal.shape[0])):
        interpolator = interpolate.interp2d(xGrid, zGrid, data_vertical[t,:,:], kind=parameterFile["interpolate"])
        #evaluate this every where to get the new interpolated grid 
        interpolatedPlaneVertical[t,:,:] = interpolator(newXGrid,newZGrid)
    #test it by plotting what the interpolated plane looks like
    #test code to plot the above and verify that it functions correctly:
    import matplotlib.pyplot as plt
    #identify the time step with the maxium average value and plot there
    t=np.argmax(data_horrizontal[:,0,0]) #time step to plot at
    plt.subplot(221)
    plt.imshow(interpolatedPlaneHorrizontal[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Horriz. Interp. '+ parameterFile["interpolate"] + " when t=" + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #now do the other one
    plt.subplot(222)
    plt.imshow(data_horrizontal[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Orig. Horriz. at time=' + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #plot the vertical plane itnerpolations
    plt.subplot(223)
    plt.imshow(interpolatedPlaneVertical[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Vertical Interp. '+ parameterFile["interpolate"] + " when t=" + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #now do the other one
    plt.subplot(224)
    plt.imshow(data_vertical[t,:,:], extent=(0,1,0,1), origin='lower')
    plt.title('Vert. orig. at time=' + str(t))
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    #save it
    plt.savefig("origionalVsInterpolated.png",dpi=400)
    plt.close()
    
    
    #assign data_horrizontal and data_vertical to the interpolated data
    data_horrizontal,data_vertical = interpolatedPlaneHorrizontal,interpolatedPlaneVertical
    
else:
    print("not interpolating!")
    target_h = h

#integrate all velocity data to form displacement data set
data_horrizontal_disp = scipy.integrate.cumtrapz(y=data_horrizontal[:,:,:],dx=dt,initial=0,axis=0)
data_vertical_disp = scipy.integrate.cumtrapz(y=data_vertical[:,:,:],dx=dt,initial=0,axis=0)
#differentiate all displacement data to form acceleration
data_horrizontal_acc = np.gradient(data_horrizontal[:,:,:],dt,axis=0)
data_vertical_acc = np.gradient(data_vertical[:,:,:],dt,axis=0)

"""

"""
Check that the integration and differentiation look reasonalbe for a point
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
ax.plot(time, data_horrizontal_disp[:,0,0])
fig.savefig("disp_plot_test.png")

fig, ax = plt.subplots()
ax.plot(time, data_horrizontal_acc[:,0,0])
fig.savefig("acc_plot_test.png")
"""
"""
#save each coordinate to the container
coordinateList = coordinates
coordinateList = list(coordinateList.reshape(coordinateList.shape[0]//3,-1,3))
#coordinateList = sorted(coordinateList,key=lambda coordinateList: coordinateList[0][2])   
for coordinate in range(len(coordinateList)): 
    #check which way is horrizontal
    if(directivity == 'sw4_x'):
        #fix the values if I have a negetive in my expression
        #checkValues = lambda x:0 if x != 0 else x     
        x0 = np.int32(np.float32(parameterFile['x0'])*precision)
        essiXstart,essiYstart =  np.int32(np.float32(parameterFile['essiXstart'])*precision),np.int32(np.float32(parameterFile['essiYstart'])*precision)
        loc_i,loc_k = int((coordinateList[coordinate][0][1]-essiYstart)/target_h)+int(x0/target_h),int((int(parameterFile['essi_z_end'])-coordinateList[coordinate][0][2])/target_h)
        #z stop if z is ever negetive!
        #assert loc_k >= 0 
        #note that j is unuesd !!
        Accelerations_dset[coordinate*3,:] = 0.0
        Accelerations_dset[coordinate*3+1,:] = data_horrizontal_acc[:,loc_i,loc_k] 
        Accelerations_dset[coordinate*3+2,:] = data_vertical_acc[:,loc_i,loc_k]  
        
        Displacements_dset[coordinate*3,:] = 0.0
        Displacements_dset[coordinate*3+1,:] = data_horrizontal_disp[:,loc_i,loc_k] 
        Displacements_dset[coordinate*3+2,:] = data_vertical_disp[:,loc_i,loc_k] 
        # print("Coordinate pair " + str(coordinateList[coordinate][0]) + " computed at " + str((loc_i,loc_k)) + "assinged to " + str(coordinate))
    
    else:
        #note that i is unuesd !!
        #checkValues = lambda x:0 if x != 0 else x      
        y0=np.int32(np.float32(parameterFile['y0'])*precision)
        essiXstart,essiYstart =  np.int32(np.float32(parameterFile['essiXstart'])*precision),np.int32(np.float32(parameterFile['essiYstart'])*precision)
        #sw4x,sw4y,sw4z = (sw4Xstart+(OSx-essixStart))/spacing,(sw4Ystart+(OSy-essiyStart))/spacing,(sw4Zstart+(essizStart-OSz))/spacing

        loc_j,loc_k = int((coordinateList[coordinate][0][0]-essiXstart)/target_h)+int(y0/target_h),int((int(parameterFile['essi_z_end'])-coordinateList[coordinate][0][2])/target_h)
        #z stop if z is ever negetive!
        #assert loc_k >= 0
        Accelerations_dset[coordinate*3,:] = data_horrizontal_acc[:,loc_j,loc_k]  
        Accelerations_dset[coordinate*3+1,:] = 0.0
        Accelerations_dset[coordinate*3+2,:] = data_vertical_acc[:,loc_j,loc_k] 
        
        Displacements_dset[coordinate*3,:] = data_horrizontal_disp[:,loc_j,loc_k]  
        Displacements_dset[coordinate*3+1,:] = 0.0
        Displacements_dset[coordinate*3+2,:] = data_vertical_disp[:,loc_j,loc_k] 
        # print("Coordinate pair " + str(coordinateList[coordinate][0]) + " computed at " + str((loc_j,loc_k)) + "assinged to " + str(coordinate))

        
#close the container
sw4essiout.close()
drm_file.close()

"""

"""
def main(args):
    #load up the parameter csv file
    
    
    
    return 0

if __name__ == '__main__':
    #main()
    pass
"""
