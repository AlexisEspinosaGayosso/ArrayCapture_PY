#!/usr/bin/env python3
#####!/usr/bin/env python

#Script for reading the different cylinder volumetric samples in time and analysing them

#Defining environment with native utilities
import os
#import math
import numpy as np
#import scipy
#import cython

FS=os.sep

#Defining environment with parallel utilities
#from joblib import Parallel, delayed
import multiprocessing

#Defining environment with my utilities
from strings_AEG import *


#defining functions that will later be integrated into another file
def unique_rows(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    return a[ui]


#The case to analyse
caseDir="/home/espinosa/ExternalMounts/IRDS/PC_ParticleCapture/disk_YY_MCPC/OpenFOAM/espinosa-2.3.x/run/FlowFields/Inline/2D/2D_inline64_SF8.0_Re100"
#caseDir="/scratch/pawsey0106/espinosa/OpenFOAM/espinosa-2.3.x/run/FlowFields/Inline/2D/2D_inline64_SF8.0_Re100"
volumetricDir=os.path.join(caseDir,"postProcessing","volumetricCirclesReal")
####For testing volumetricDir=os.path.join(caseDir,"postProcessing","test")

#The collectors centres in the basic mesh
lowLeft=np.array([-0.015666,  -0.015666,  -8.8620278672547E-4])
topRight=np.array([0.046998,  0.046998,  8.8620278672547E-4])
deltaL=topRight-lowLeft
DCollector=0.01
rCollector=DCollector/2
cylNamesBase=list(["cylinder_downLeft","cylinder_upLeft", "cylinder_upRight", "cylinder_downRight" ])
cylCentresBase=np.array([[0, 0], [0, 0.031332], [0.031332, 0.031332], [0.031332, 0]])
NCylBase=len(cylNamesBase)

#The python processed names
pythonDir=os.path.join(caseDir,"postProcessing","pythonFiles")
if not os.path.exists(pythonDir):
    os.makedirs(pythonDir)
outVolDir=os.path.join(pythonDir, "Volumetric")
if not os.path.exists(outVolDir):
    os.makedirs(outVolDir)

#The existing times in the volumetric directory
listTimesRaw=os.listdir(volumetricDir)
#print(listTimesRaw)
listTimesClean=sorted( list( [ x for x in listTimesRaw if is_number(x)]))
#print(listTimesClean)
arrayTimes=np.array(list(map(float, listTimesClean)))
Ntimes=len(arrayTimes)
arrayTimesSorted=sorted(list(arrayTimes))
if np.array_equal(arrayTimesSorted, arrayTimes):
    print("Time array is ready and has",  Ntimes,  " times")
else:
    print("Time array from directories has a sorting problem")

#The cylinder names array
#Intended to read the cylinder names directly from the files:
####listCylNamesRaw=os.listdir(volumetricDir+"/"+listTimesClean[0])
####listCylNamesClean=list([x for x in listCylNamesRaw if "cylinder" in x ])
####Ncylinders=len(listCylNamesClean)
####print("There are ", Ncylinders, " cylinders")
#Generating the list of cylinder names and the array of their centres with a cycle
NSX=4
NSY=4
####NSX=1 #For testing
####NSY=1
listCylNamesClean=list([]);
centresClean=np.reshape([-1000, -1000], (1, 2)); #dummy trick to start with an array of the right shape
print(centresClean.shape)
for ii in range(0, NSX):
    for jj in range(0, NSY):
        for kk in range(0, NCylBase):
            listCylNamesClean.append("rp009008_"+cylNamesBase[kk]+"_piece_"+str(ii+1)+"_"+str(jj+1)+"_Real_U.xy") 
            #centresClean=np.append(centresClean, cylCentresBase[kk, :]+np.array([ii*deltaL[0], jj*deltaL[1]]), axis=0)
            newCentre=np.reshape(cylCentresBase[kk, :], (1, 2))+np.reshape([ii*deltaL[0], jj*deltaL[1]], (1, 2))
            centresClean=np.append(centresClean, newCentre, axis=0)
centresClean=np.delete(centresClean, 0, 0)
print(centresClean.shape)
print(centresClean)
NCylinders=len(centresClean[:, 1])
print("There are NCylinders=",  NCylinders)



#------------------------------
#Cycle for processing each cylinder


#New for parallel. Preparation of the input for the new function calcOneVolumetric,
inputs = [];
for ii in range(0,NCylinders):
    inputs.append((ii))
    
#Two New lines for the definition a function that will be operated in parallel instead of the original serial for loop,
def calcOneVolumetric(iiAndOthers):
    ii = iiAndOthers

#Original line defining the serial for loop,
####for ii in range(0,NCylinders):

#Rest of the lines are the same, either for the for loop or the function (defined above)
    iName=listCylNamesClean[ii]
    outputName,  extension =os.path.splitext( iName)
    outputName=os.path.join(outVolDir, outputName+"_VolumetricFlux"+".dat") 
    print(outputName)
    existingTime=np.reshape([-1000], (1, 1)); #dummy trick to start with an array of the right shape
    if os.path.isfile(outputName):
        MAux=np.genfromtxt(outputName)
        if len(MAux)>0:
            existingTime=np.genfromtxt(outputName)[:, 0]
        #End if len
    #End if os.path.isfile
    existingLast=existingTime[-1]
    for jj in range(0,Ntimes):        
        jTimeC=listTimesClean[jj]
        jTimeN=arrayTimes[jj]
        if jTimeN > existingLast:
            dataName=os.path.join(volumetricDir, jTimeC, iName)        
            #print(dataName)
            dataArray=np.genfromtxt(dataName, delimiter=None) #, usecols=(0, 1)) 
            
            #Reading the coordinates and values, and converting to local polar coordinates
            xLocal=dataArray[:, 0]-centresClean[ii, 0]
            yLocal=dataArray[:, 1]-centresClean[ii, 1]
            zLocal=dataArray[:, 2]
            thetaLocal=np.arctan2(yLocal, xLocal)
            rpLocal=(np.sqrt(np.power(xLocal, 2)+np.power(yLocal, 2))-rCollector)/rCollector
            UxLocal=dataArray[:, 3]
            UyLocal=dataArray[:, 4]
            UzLocal=dataArray[:, 5]
            UrLocal=np.multiply(UxLocal, np.cos(thetaLocal))+np.multiply(UyLocal, np.sin(thetaLocal))
            UthetaLocal=-np.multiply(UxLocal, np.sin(thetaLocal))+np.multiply(UyLocal, np.cos(thetaLocal))
                      
            #Sorting the arrays with the thetaLocal
            zipA=np.stack((thetaLocal, rpLocal, xLocal, yLocal, zLocal, UrLocal, UthetaLocal), 1)
            zipSort=zipA[np.argsort(zipA[:, 0])]
            deltaTheta=zipSort[1:, 0]-zipSort[0:-1, 0]
            deltaTheta=np.reshape(np.insert(deltaTheta, 0, np.amax(deltaTheta)), (len(deltaTheta)+1, 1))
            zipSort=np.append(zipSort, deltaTheta,axis=1)        
            lineOut=np.argmin(zipSort, axis=0)[-1]
            zipSort=np.delete(zipSort, lineOut, axis=0)
            
            
            #Integrating the fluxes towards the cylinder
            UrLocal=zipSort[:, -3] #Third column from end to beginning is Ur
            negativeValues=UrLocal[UrLocal<0]
            kArea=deltaTheta[0]*(rpLocal[0]+1)*rCollector*deltaL[2]
            volumetricFlux=-np.sum(negativeValues)*kArea
            outputLine=np.reshape([jTimeN, volumetricFlux], (1, 2))

            #Writing the new volumetric integral into this time
            fOut=open(outputName, mode='ab')
            np.savetxt(fOut, outputLine)
            fOut.close()             
        #End if jTimeN not in
    #End for jj in range    
    prueba=np.genfromtxt(outputName)
    print("Aqui tienes la prueba")
    print(prueba) 
#End of the for ii in range(0,NCylinders) .OR. def calcOneVolumetric function

#Declaring the number of parallel workers
num_of_workers = multiprocessing.cpu_count()
pool = multiprocessing.Pool(num_of_workers)

#Operating the calcOneVolumetric function in parallel
for data in pool.map(calcOneVolumetric, inputs):
    ii=data
        
#End of script
print("Script done")
    

