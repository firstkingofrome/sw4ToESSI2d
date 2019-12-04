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

def roundToBase(number,base):
    return np.int32(base*np.round(number/base))


#def main(args):
"""
sw4essiout_filename =  args.sw4HDF5Out[0]
pointsRequested = args.refPoints[0]
"""
sw4essiout_filename="Location17faultParallel.cycle=00000.essi"
pointsRequested="Comparison_points.csv"
pointsOutput = "outputPoints"
#load the csv file  
points = pandas.read_csv(pointsRequested)
#open the hdf5 container for this run and get each out put velocity point
sw4essiout = h5py.File(sw4essiout_filename, 'r')
#get constant values
sw4Xstart,sw4Ystart,sw4Zstart = int(points["sw4Xstart"][0]),int(points["sw4Ystart"][0]),int(points["sw4Zstart"][0])
essixStart,essiyStart,essizStart = int(points["essiXstart"][0]),int(points["essiYstart"][0]),int(points["essiZstart"][0])
spacing = int(sw4essiout['ESSI xyz grid spacing'][:][0])
print("rounding all data to nearest " + str(spacing) + "'s")
#fix the precision of everything based on the spacing specified by the sw4 model
base = spacing
essixStart,essiyStart,essizStart = roundToBase(essixStart,base),roundToBase(essiyStart,base),roundToBase(essizStart,base)

dt = sw4essiout['timestep'][0]
"""
Used for plot unit tests
"""
npts = sw4essiout['vel_0 ijk layout'].shape[0]
t1 = dt*(npts-1)
t0 = sw4essiout['time start'][0]
time = np.linspace(t0, t1, npts)
print("CHECK THE SPACING SW4 CURRENTLY THINKS IT IS " +str(spacing))
#check directivity
if(points["azimuth"][0]==90):
    #I am using the sw4 x plane direction |
    directivity = 'sw4_x'
else:
    #I am using the sw4 y plane direction --
    directivity = 'sw4_y'

#spacing = int(points["Spacing"][0])
for index, row in points.iterrows():
    print(row)
    #translate into xyz coordinates
    ### THE ESSI CORDINATES ARE XYZ WHILE SW4 IS YXZ
    OSx,OSy,OSz = roundToBase(row["x"],base),roundToBase(row["y"],base),roundToBase(row["z"],base)
    #round all of this to the nearest point as defined by the spacing
    print("preparing refrence point " + row["point"])
    print("saving x y and z velocities for  exssi point " + "%i,%i,%i" % (OSx,OSy,OSz))
    """
    THE LINE BELOW THIS IS WHERE I CONVERT FROM ESSI --> SW4
    """
    if(directivity == 'sw4_x'):
        sw4x,sw4y,sw4z = (sw4Xstart+(OSy-essiyStart))/spacing,(sw4Ystart+(OSx-essixStart))/spacing,(sw4Zstart+(essizStart-OSz))/spacing
        print("This corresponds to sw4 index %i,%i,%i"  %(sw4x,sw4y,sw4z))
        print('\n')
        #save xyz as csv files for that point
        #save ESSI x y and z!
        fname = "_sw4_"+row["point"]+".csv"
        x = sw4essiout['vel_0 ijk layout'][:,sw4x,sw4y,sw4z]
        np.savetxt(fname='x'+"_vel"+fname,X=x)
        #save y
        #y = sw4essiout['vel_0 ijk layout'][:,sw4x,sw4y,sw4z]
        #np.savetxt(fname='y'+fname,X=y)
        #save z
        z = -1.0*sw4essiout['vel_2 ijk layout'][:,sw4x,sw4y,sw4z]
        np.savetxt(fname='z'+"_vel"+fname,X=z)

    else:
        sw4x,sw4y,sw4z = (sw4Xstart+(OSy-essiyStart))/spacing,(sw4Ystart+(OSx-essixStart))/spacing,(sw4Zstart+(essizStart-OSz))/spacing
        print("This corresponds to sw4 index %i,%i,%i"  %(sw4x,sw4y,sw4z))
        print('\n')
        #save xyz as csv files for that point
        #save ESSI x y and z!
        fname = "_sw4_"+row["point"]+".csv"
        x = sw4essiout['vel_1 ijk layout'][:,sw4x,sw4y,sw4z]
        np.savetxt(fname='x'+"_vel"+fname,X=x)
        #save y
        #y = sw4essiout['vel_0 ijk layout'][:,sw4x,sw4y,sw4z]
        #np.savetxt(fname='y'+fname,X=y)
        #save z
        z = -1.0*sw4essiout['vel_2 ijk layout'][:,sw4x,sw4y,sw4z]
        np.savetxt(fname='z'+"_vel"+fname,X=z)



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
