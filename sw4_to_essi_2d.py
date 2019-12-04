#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  get2dSliceNP.py
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

#load the input data
essiXstart,essiYStart,essiZstart=int(parameterFile['essiXstart']),int(parameterFile['essiYstart']),int(parameterFile['essiZstart'])
plane = int(int(parameterFile['horizontal_plane'])/h)

if(directivity == 'sw4_x'):
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



#load up the geometery file
drm_file = h5py.File(drm_filename,"a")
coordinates = drm_file['Coordinates']
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

#integrate all velocity data to form displacement data set
data_horrizontal_disp = scipy.integrate.cumtrapz(y=data_horrizontal[:,:,:],dx=dt,initial=0,axis=0)
data_vertical_disp = scipy.integrate.cumtrapz(y=data_vertical[:,:,:],dx=dt,initial=0,axis=0)
#differentiate all displacement data to form acceleration
data_horrizontal_acc = np.gradient(data_horrizontal[:,:,:],dt,axis=0)
data_vertical_acc = np.gradient(data_vertical[:,:,:],dt,axis=0)

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
#save each coordinate to the container
coordinateList = coordinates[:]
coordinateList = list(coordinateList.reshape(coordinateList.shape[0]//3,-1,3))
#coordinateList = sorted(coordinateList,key=lambda coordinateList: coordinateList[0][2])   
for coordinate in range(len(coordinateList)): 
    #check which way is horrizontal
    if(directivity == 'sw4_x'):
        #fix the values if I have a negetive in my expression
        #checkValues = lambda x:0 if x != 0 else x         
        loc_i,loc_j,loc_k = int(coordinateList[coordinate][0][0]/h),int(coordinateList[coordinate][0][1]/h),int((essiZstart-coordinateList[coordinate][0][2])/h)
        #z stop if z is ever negetive!
        #assert loc_k >= 0 
        #note that j is unuesd !!
        Accelerations_dset[coordinate*3,:] =data_horrizontal_acc[:,loc_i,loc_k] 
        Accelerations_dset[coordinate*3+1,:] = 0.0
        Accelerations_dset[coordinate*3+2,:] = data_vertical_acc[:,loc_i,loc_k]  
        
        Displacements_dset[coordinate*3,:] = data_horrizontal_disp[:,loc_i,loc_k] 
        Displacements_dset[coordinate*3+1,:] = 0.0
        Displacements_dset[coordinate*3+2,:] = data_vertical_disp[:,loc_i,loc_k] 
        print("Coordinate pair " + str(coordinateList[coordinate][0]) + " computed at " + str((loc_i,loc_k)) + "assinged to " + str(coordinate))
    
    else:
        #note that i is unuesd !!
        #checkValues = lambda x:0 if x != 0 else x         
        loc_j,loc_k = int((coordinateList[coordinate][0][0]-essiXstart)/h),int((int(parameterFile['essi_z_end'])-coordinateList[coordinate][0][2])/h)
        #z stop if z is ever negetive!
        #assert loc_k >= 0
        Accelerations_dset[coordinate*3,:] = data_horrizontal_acc[:,loc_j,loc_k]  
        Accelerations_dset[coordinate*3+1,:] = 0.0
        Accelerations_dset[coordinate*3+2,:] = data_vertical_acc[:,loc_j,loc_k] 
        
        Displacements_dset[coordinate*3,:] = data_horrizontal_disp[:,loc_j,loc_k]  
        Displacements_dset[coordinate*3+1,:] = 0.0
        Displacements_dset[coordinate*3+2,:] = data_vertical_disp[:,loc_j,loc_k] 
        print("Coordinate pair " + str(coordinateList[coordinate][0]) + " computed at " + str((loc_j,loc_k)) + "assinged to " + str(coordinate))

        
#close the container
#sw4essiout.close()
#drm_file.close()


"""
def main(args):
    #load up the parameter csv file
    
    
    
    return 0

if __name__ == '__main__':
    #main()
    pass
"""
