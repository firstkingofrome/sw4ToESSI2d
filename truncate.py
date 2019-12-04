#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import h5py
import os
import numpy as np
inputHDF5File = "DRM_input_Model27.hdf5"
outputHDF5File = "TRUNCATED_DRM_input_Model27.hdf5"
startTime,stopTime=6,20
#make a copy of the inputHDF5File
#os.system("rm " + outputHDF5File)
os.system("cp "+inputHDF5File +" " + outputHDF5File)
#load up the file
fileToTruncate = h5py.File(outputHDF5File,"r+")
#compute the time range to truncate over
dt = fileToTruncate["Time"][1]-fileToTruncate["Time"][0]
startTimeStep,stopTimeStep = int(startTime/dt),int(stopTime/dt)
#truncate time
time = fileToTruncate["Time"][startTimeStep:stopTimeStep]
#compute new time range starting at 0
time = np.linspace(0,time[-1]-startTime,time.shape[0])
del fileToTruncate["Time"]
fileToTruncate.create_dataset("Time",data=time)
#now truncate accelerations and displacements
Accelerations = fileToTruncate["Accelerations"][:,startTimeStep:stopTimeStep]
Displacements = fileToTruncate["Displacements"][:,startTimeStep:stopTimeStep]
del fileToTruncate["Displacements"]
del fileToTruncate["Accelerations"]
fileToTruncate.create_dataset("Accelerations",data=Accelerations)
fileToTruncate.create_dataset("Displacements",data=Displacements)
