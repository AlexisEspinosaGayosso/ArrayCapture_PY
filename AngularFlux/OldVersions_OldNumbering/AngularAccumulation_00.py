#!/usr/bin/env python3
#####!/usr/bin/env python

#Script for reading the different cylinder volumetric samples in time and  performing angular accumulation
#Defining environment with native utilities
import os
import shutil
import sys
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

#defining a function to delete the last line of a file
def deleteLastLine(fileName):
    afile=open(fileName,"rb+") #Opening the file in binary mode
    afile.seek(0, os.SEEK_END) #Move the pointer to the end of the file
    pos = afile.tell() -1 #Define the position one character before the end
    while pos > 0 and afile.read(1) != b"\n": #Read backwards until finding another EOL
        pos -= 1                                                   #Note that b"\n" is the binary character EOL
        afile.seek(pos,  os.SEEK_SET)
    #End while pos >0
    if pos > 0: #Truncate everything from where we are using truncate
        afile.seek(pos,  os.SEEK_SET)
        afile.truncate()
        afile.write(b"\n") #Adding an end of line after truncating because I need this character at the end of the exisiting line
    #End if pos>0    
    afile.close()
#End def deleteLastLine


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
outAngularDir=os.path.join(pythonDir, "Angular")
if not os.path.exists(outAngularDir):
    os.makedirs(outAngularDir)

#The existing times in the volumetric directory
listTimesRaw=os.listdir(volumetricDir)
print(listTimesRaw[-1])
##listTimesClean=sorted( list( [ x for x in listTimesRaw if is_number(x)]))
listTimesClean=list( [ x for x in listTimesRaw if is_number(x)])
#print(listTimesClean[-1])
listTimesClean.sort(key=float)
#print(listTimesClean[-1])
arrayTimes=np.array(listTimesClean, dtype=np.float)
#print(arrayTimes[-1])
#print(np.shape(arrayTimes))
Ntimes=len(arrayTimes)
arrayTimesSorted=np.sort(arrayTimes)
#print(arrayTimesSorted[-1])
#print(np.shape(arrayTimesSorted))
if np.array_equal(arrayTimesSorted, arrayTimes):
    print("Time array is ready and has",  Ntimes,  " times")
else:
    print("Time array from directories has a sorting problem")
    listFile=open(os.path.join(outAngularDir, "aaatimesSortedText.dat"),'w')
    timesFile=open(os.path.join(outAngularDir, "aaatimesSortedNumerical.dat"),'w')
    np.savetxt(listFile, arrayTimes)
    np.savetxt(timesFile, arrayTimesSorted)
    listFile.close()  
    timesFile.close()
    sys.exit()
#End if np.array_equal

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
    baseName,  extension =os.path.splitext( iName)
    outputName=os.path.join(outAngularDir, baseName+"_AngularFlux"+".dat") 
    backupName=os.path.join(outAngularDir, baseName+"_AngularFlux"+".dat"+".bak") 
    print(outputName)
    existingTime=np.reshape([-1000, -1000], (1, 2)); #dummy trick to start with an array of the right shape
    if os.path.isfile(outputName):
        try:
            print("Reading the matrix for the first time, cylinder: ",  baseName)
            MAux=np.genfromtxt(outputName)
        except:
            print("Something went wrong with first reading, cylinder: ",  baseName)
            print("Backing up to file .back0 and Deleting the last line, cylinder: ",  baseName)
            shutil.copy(outputName, backupName+"0")
            deleteLastLine(outputName)
            print("Reading the matrix again, cylinder: ",  baseName)
            MAux=np.genfromtxt(outputName)
        #End try
        if len(MAux)>0:
            try:
                print("Reading existingTime, cylinder: ",  baseName)
                existingTime=np.genfromtxt(outputName)[:, [0, 1]]
            except:
                print("Something went wrong with existingTime reading, cylinder: ",  baseName)
                print("Backing up to file .back1 and Deleting the last line, cylinder: ",  baseName)
                shutil.copy(outputName, backupName+"1")
                deleteLastLine(outputName)
                print("Reading the existingTime again, cylinder: ",  baseName)
                existingTime=np.genfromtxt(outputName)[:, [0, 1]]
            #End try
        #End if len
    #End if os.path.isfile
    existingLastEnd=existingTime[-1, 1]
    existingLastIni=existingTime[-1, 0]
    arrayTimeLast=arrayTimes[-1]    
    print("existingLast: ",  [existingLastIni, existingLastEnd])
    if (existingLastEnd<arrayTimeLast):
        startIndexTimes=int(np.searchsorted(arrayTimes,existingLastIni)) #Starting from the  beginning of the last line
        if os.path.isfile(outputName):
            shutil.copy(outputName, backupName)
            deleteLastLine(outputName) #Deleting the last line in case it is corrupted
            print("Before starting posprocess, Deleting last line of the result file: ",  baseName)
        #End if os.path.isfile
        
        #------------------------------------------------------------------------------------------------------------------------
        #Initial reading of first data and ordering of angular samples ------------------------------
        jTimeCS=listTimesClean[startIndexTimes]
        jTimeNS=arrayTimes[startIndexTimes]
        print("Starting postprocess on time: ",  jTimeNS)
        dataNameS=os.path.join(volumetricDir, jTimeCS, iName)
        dataArrayS=np.genfromtxt(dataNameS, delimiter=None) #, usecols=(0, 1)) 
                    
        #Reading the coordinates and values, and converting to local polar coordinates
        xLocalS=dataArrayS[:, 0]-centresClean[ii, 0]
        yLocalS=dataArrayS[:, 1]-centresClean[ii, 1]
        zLocalS=dataArrayS[:, 2]
        thetaLocalS=np.arctan2(yLocalS, xLocalS)
        rpLocalS=(np.sqrt(np.power(xLocalS, 2)+np.power(yLocalS, 2))-rCollector)/rCollector
        UxLocalS=dataArrayS[:, 3]
        UyLocalS=dataArrayS[:, 4]
        UzLocalS=dataArrayS[:, 5]
        UrLocalS=np.multiply(UxLocalS, np.cos(thetaLocalS))+np.multiply(UyLocalS, np.sin(thetaLocalS))
        UthetaLocalS=-np.multiply(UxLocalS, np.sin(thetaLocalS))+np.multiply(UyLocalS, np.cos(thetaLocalS))
                  
        #Sorting the arrays with the thetaLocal
        zipAS=np.stack((thetaLocalS, rpLocalS, xLocalS, yLocalS, zLocalS, UrLocalS, UthetaLocalS), 1)
        sortingOrderS=np.argsort(zipAS[:, 0])
        zipSortS=zipAS[sortingOrderS]
        deltaThetaS=zipSortS[1:, 0]-zipSortS[0:-1, 0]
        deltaThetaS=np.reshape(np.insert(deltaThetaS, 0, np.amax(deltaThetaS)), (len(deltaThetaS)+1, 1))
        zipSortPlusS=np.append(zipSortS, deltaThetaS,axis=1)
        lineOutS=np.argmin(zipSortPlusS, axis=0)[-1]
        zipSortMinusS=np.delete(zipSortPlusS, lineOutS, axis=0)
          
        #Integrating the fluxes towards the cylinder
        UrLocalS=zipSortMinusS[:, -3] #Third column from end to beginning is Ur
        rpLocalS=zipSortMinusS[:,1] #Second column (starting from 0 index) is rpLocal
        thetaLocalS=zipSortMinusS[:,0] #First column (starting from 0 index) is thetaLocal
        NTheta=np.size(thetaLocalS)
        UrMClippedS=np.negative(np.clip(UrLocalS, -np.finfo('d').max, 0))
        ##kAreaS=deltaThetaS[0]*(rpLocalS[0]+1)*rCollector*deltaL[2]
        ##radialVolumeAccS=UrMClippedS*kAreaS*deltaT    
        #------------------------------------------------------------------------------------------------------------------------

        #Writing the angles headers if it is the first line to write
        if (existingLastEnd == -1000):
            outputLine=np.reshape(np.insert(thetaLocalS,0,  [-1000, -1000]), (1, NTheta+2))
            fOut=open(outputName, mode='ab')
            np.savetxt(fOut, outputLine)
            fOut.close()
              

        #Finding the beginning and end of accumulation for the first time
        indexIni=startIndexTimes
        timeIni=arrayTimes[indexIni]
        timeEnd=np.floor(timeIni+1) #Tentative
        indexEnd=int(np.searchsorted(arrayTimes,timeEnd)) 
        indexEnd=np.minimum(indexEnd, Ntimes-1)
        timeEnd=arrayTimes[indexEnd]
        
        #Declaring the array for accumulation starting from zero
        radialVolumeAcc=np.zeros(UrMClippedS.shape)
        
        for jj in range(startIndexTimes+1,Ntimes):           
            jTimeC=listTimesClean[jj]
            jTimeN=arrayTimes[jj]
            deltaT=jTimeN-jTimeNS
            
            
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
            sortingOrder=np.argsort(zipA[:, 0])
            zipSort=zipA[sortingOrder]
            deltaTheta=zipSort[1:, 0]-zipSort[0:-1, 0]
            deltaTheta=np.reshape(np.insert(deltaTheta, 0, np.amax(deltaTheta)), (len(deltaTheta)+1, 1))
            zipSortPlus=np.append(zipSort, deltaTheta,axis=1)        
            lineOut=np.argmin(zipSortPlus, axis=0)[-1]
            zipSortMinus=np.delete(zipSortPlus, lineOut, axis=0)
                        
            #Integrating the fluxes towards the cylinder
            UrLocal=zipSortMinus[:, -3] #Third column from end to beginning is Ur
            rpLocal=zipSortMinus[:,1] #Second column (starting from 0 index) is rpLocal
            UrMClipped=np.negative(np.clip(UrLocal, -np.finfo('d').max, 0))
            kArea=deltaTheta[0]*(rpLocal[0]+1)*rCollector*deltaL[2]
            radialVolumeAcc=radialVolumeAcc+(UrMClipped+UrMClippedS)/2.0*kArea*deltaT        

            #Final operations in the cycle
            
            if ((jj == indexEnd) or (jj == Ntimes-1)):
                #Writing the new volumetric integral into this time
                outputLine=np.reshape(np.insert(radialVolumeAcc,0,  [timeIni, jTimeN]), (1, NTheta+2))
                fOut=open(outputName, mode='ab')
                np.savetxt(fOut, outputLine)
                fOut.close()
            #End if jj==indexEnd
            
            if ((jj == indexEnd) and (jj != Ntimes-1)):    
                #Finding the beginning and end of accumulation for the first time
                indexIni=indexEnd #Starting from the  beginning of the last line
                timeIni=arrayTimes[indexIni]
                timeEnd=np.floor(timeIni+1) #Tentative
                indexEnd=int(np.searchsorted(arrayTimes,timeEnd))
                indexEnd=np.minimum(indexEnd, Ntimes-1)
                timeEnd=arrayTimes[indexEnd]
                
                #Declaring the array for accumulation starting from zero
                radialVolumeAcc=np.zeros(UrMClippedS.shape)
            #End if jj==indexEnd
            
            #Cyclic substitution
            UrMClippedS=UrMClipped
            jTimeNS=jTimeN
        #End for jj in range
    #End if  existingLastEnd<arrayTimeLast  
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
    

