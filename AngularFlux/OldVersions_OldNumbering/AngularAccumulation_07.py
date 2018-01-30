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
import csv

FS=os.sep

#Defining print precision
np.set_printoptions(precision=15)

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

#defining a function to copy a file catching error messages
def copyFile(src, dest):
    try:
        shutil.copy(src, dest)
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)
#End def copyFile        


#The case to analyse
#caseDir="/home/espinosa/ExternalMounts/IRDS/PC_ParticleCapture/disk_YY_MCPC/OpenFOAM/espinosa-2.3.x/run/FlowFields/Inline/2D/2D_inline64_SF8.0_Re100"
#caseDir="/home/espinosa/ExternalMounts/IRDS/PC_ParticleCapture/disk_YY_MCPC/OpenFOAM/espinosa-2.3.x/run/FlowFields/Inline/2D/2DFull_inline120Alexis_SF8.0_Re100"
#caseDir="/scratch/pawsey0106/espinosa/OpenFOAM/espinosa-2.3.x/run/FlowFields/Inline/2D/2D_inline64_SF8.0_Re100"
#volumetricDir=os.path.join(caseDir,"postProcessing","volumetricCirclesReal")
####For testing volumetricDir=os.path.join(caseDir,"postProcessing","test")

#Format for the rp in file names
formatRp="%10.8f"

#Checking the arguments
NArgNeeded=4
if len(sys.argv) != NArgNeeded+1:
    print("This script needs " + str(NArgNeeded) + " arguments.")
    print("The first is the rp size")
    print("The second is the Solid Fraction SF")
    print("The third is the Mean or Normal variables")
    print("The fourth is the directory caseDir")
    print("So a usage example goes like this:")
    print("python " + sys.argv[0] + " 0.0053091 8.0 ./2DFull_inline")
    print("Exiting the script")
    sys.exit(1)
else:
    try:
        #The rp to be analysed
        rpAnalysis=np.float64(sys.argv[1])
    except:
        print("The rp size is not a float number")
        print("You gave: ",  sys.argv[1])
        print("Exiting the script")
        sys.exit(1)
    try:
        #The SF to be used
        SFAnalysis=np.float64(sys.argv[2])
    except:
        print("The SFAnalysis is not a float number")
        print("You gave: ",  sys.argv[2])
        print("Exiting the script")
        sys.exit(1)
    #The caseDir to be used
    caseDir=sys.argv[4]
    if not os.path.exists(caseDir):
        print("The caseDir directory does not exists")
        print("You gave: ",  sys.argv[4])
        print("Exiting the script")
        sys.exit(1)
    else:
        #The directory to read
        tagMean=""
        if (sys.argv[3]=="Mean" or sys.argv[3]=="mean"):
            tagMean="Mean"
        volumetricDir=os.path.join(caseDir,"postProcessing","velocity"+tagMean+"Rings")
        if not os.path.exists(volumetricDir):
            print("The directory")
            print(volumetricDir)
            print(" does not exists")
            print("You gave: ",  sys.argv[3],  "as third argument")
            print("Exiting the script")
            sys.exit(1)
        #End if not os.path.exists(volumetricDir)
    #End if not os.path.exisits(caseDir)
#End if len(sys.argv)
                
#Reading the content of the geometryFile with csv
geometryFile=os.path.join(caseDir,"postProcessing", "generalGeometry.txt");
if not os.path.isfile(geometryFile):
    print("The geometry file does not exist")
    print("We are looking for: ", geometryFile)
    print("Exiting the script")
    sys.exit(2);
else:
    ifile = open(geometryFile, "r")
    reader=csv.reader(ifile)
    for row in reader:
        if row[0] == "ArrayType":
            arrayType=row[1]
            print(arrayType)
        elif row[0] == "NCylBlock":
            NCylBlock=np.int(row[1])
            print(NCylBlock)
        elif row[0] == "NBX":
            NSX=np.int(row[1])
            print(NSX)
        elif row[0] == "NBY":
            NSY=np.int(row[1])
            print(NSY)
        elif row[0] == "SF":
            SFArray=np.array(row[1:], dtype=np.float64)
            print(SFArray)
        elif row[0] == "LOneCyl":
            LOneCylArray=np.array(row[1:], dtype=np.float64)
            print(LOneCylArray)
        elif row[0] == "ZOneCyl":
            ZOneCyl=np.float64(row[1])
            print(ZOneCyl)
        elif row[0] == "DC":
            DCollector=np.float64(row[1])
            rCollector=DCollector/2.0
            print(DCollector)
        elif row[0] == "NthetaDiv":
            NthetaDiv=np.int(row[1])
            print(NthetaDiv)
        elif row[0] == "DomainType":
            domainType=row[1]
            print(domainType)
        #End the case
    #End for row
#End if not geometry file

#Correcting numbers when the simulation was a Channel
try:
    if domainType=="Channel":        
        NSY=1 #Correcting the number of rows
        NthetaDiv=np.int(np.round(NthetaDiv/2)) #Correcting the number of divisions
        print("Changing variables because of the channel process")
        print("Now NSY=",  NSY,  " and NthetaDiv=", NthetaDiv)
    #End if domainType==Channel
except:
    domainType=""
#End try for domainType

#The rSampling array
samplingRpFile=os.path.join(caseDir,"postProcessing", "samplingRp.txt");
if not os.path.isfile(samplingRpFile):
    print("The samplingRp file does not exist")
    print("We are looking for: ", samplingRpFile)
    print("Exiting the script")
    sys.exit(3);
else:
    rpSampling=np.genfromtxt(samplingRpFile)
    samplingI=np.int(np.searchsorted(rpSampling,rpAnalysis))
    rpHere=rpSampling[samplingI]
    rowIni=np.int(NthetaDiv*samplingI )
    rowEnd=np.int(NthetaDiv*(samplingI+1))
    print("rowIni=",  rowIni,  " , rowEnd=",  rowEnd)
#End if samplingRpFile

#The case basic length LOneCyl
iSF=np.int(np.searchsorted(SFArray,SFAnalysis))
LOneCyl=LOneCylArray[iSF]

#The collectors centres in the basic mesh
if arrayType=="Inline":
    lowLeft=np.array([-LOneCyl/2.0,  -LOneCyl/2.0,  -ZOneCyl], dtype=np.float64)
    topRight=np.array([(LOneCyl/2.0)*3.0,  (LOneCyl/2.0)*3.0,  ZOneCyl], dtype=np.float64)
    cylCentresBase=np.array([[0, 0], [0, LOneCyl], [LOneCyl, LOneCyl], [LOneCyl, 0]], dtype=np.float64)    
    cylNamesBase=list(["cylinder_downLeft","cylinder_upLeft", "cylinder_upRight", "cylinder_downRight" ])
elif arrayType=="Staggered":
    lowLeft=np.array([-(LOneCyl/2.0)*5.0,  -LOneCyl/2.0,  -ZOneCyl], dtype=np.float64)
    topRight=np.array([(LOneCyl/2.0)*3.0,  (LOneCyl/2.0)*3.0,  ZOneCyl], dtype=np.float64)
    cylCentresBase=np.array([[-LOneCyl*2.0,  0], [-LOneCyl, LOneCyl], [LOneCyl, LOneCyl],  [0, 0]], dtype=np.float64)
    cylNamesBase=list(["cylinder_downLeft","cylinder_upLeft", "cylinder_upRight", "cylinder_downRight" ])
elif arrayType=="Single":
    lowLeft=np.array([-LOneCyl,  -LOneCyl,  -ZOneCyl], dtype=np.float64)
    topRight=np.array([LOneCyl,  LOneCyl,  ZOneCyl], dtype=np.float64)
    cylCentresBase=np.array([[0, 0]], dtype=np.float64)
    cylNamesBase=list(["cylinder" ])
#End if arrayType
deltaL=topRight-lowLeft
NCylBase=len(cylNamesBase)

#The python processed names
pythonDir=os.path.join(caseDir,"postProcessing","pythonFiles")
if not os.path.exists(pythonDir):
    os.makedirs(pythonDir)
outAngularDir=os.path.join(pythonDir, "Angular"+tagMean)
if not os.path.exists(outAngularDir):
    os.makedirs(outAngularDir)

#The existing times in the volumetric directory
listTimesRaw=os.listdir(volumetricDir)
print(listTimesRaw[-1])
listTimesClean=list( [ x for x in listTimesRaw if is_number(x)])
listTimesClean.sort(key=float)
arrayTimes=np.array(listTimesClean, dtype=np.float)
Ntimes=len(arrayTimes)
arrayTimesSorted=np.sort(arrayTimes)

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
    sys.exit(4)
#End if np.array_equal

#The cylinder names array
listCylNamesClean=list([]);
centresClean=np.reshape([-1000, -1000], (1, 2)); #dummy trick to start with an array of the right shape
print(centresClean.shape)
for ii in range(0, NSX):
    for jj in range(0, NSY):
        for kk in range(0, NCylBase):
            listCylNamesClean.append(cylNamesBase[kk]+"_piece_"+str(ii+1)+"_"+str(jj+1)+"_velocity"+tagMean+"Rings_U.xy") 
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

#New for parallel. Preparation of the input for the new function calcOneAngular,
inputs = [];
for ii in range(0,NCylinders):
    inputs.append((ii))
    
#Two New lines for the definition a function that will be operated in parallel instead of the original serial for loop,
def calcOneAngular(iiAndOthers):
    ii = iiAndOthers

#Original line defining the serial for loop,
####for ii in range(0,NCylinders):

    #Defining the file and directories names
    iName=listCylNamesClean[ii]
    baseName,  extension =os.path.splitext( iName)
    outRpDir=os.path.join(outAngularDir, "Rp_" + (formatRp % rpHere))
    if not os.path.exists(outRpDir):
        try:
           os.makedirs(outRpDir)
        except:
           if not os.path.exists(outRpDir):
               print("There is still a problem in the creation of the directory:")
               print(outRpDir)
               sys.exit(5)
           #End if not exists
        #End try
    #End if not exists
              
    outputName=os.path.join(outRpDir, (formatRp % rpHere) + "_" + baseName + "_Angular"+tagMean+"Flux"+".dat") 
    backupName=os.path.join(outRpDir, (formatRp % rpHere) + "_" + baseName + "_Angular"+tagMean+"Flux" + ".dat"+".bak") 
    print(outputName)
    
    #Reading the existing file and defininig the initial time
    existingTime=np.reshape([-1000, -1000], (1, 2)); #dummy trick to start with an array of the right shape
    if os.path.isfile(outputName):
        try:
            print("Reading the matrix for the first time, cylinder: ",  baseName)
            MAux=np.genfromtxt(outputName)
        except:
            print("Something went wrong with first reading, cylinder: ",  baseName)
            print("Backing up to file .bak0 and Deleting the last line, cylinder: ",  baseName)
            copyFile(outputName, backupName+"0")
            deleteLastLine(outputName)            
            try:
                print("Reading the matrix again, cylinder: ",  baseName)
                MAux=np.genfromtxt(outputName)
            except:
                print("Something went double wrong with first reading, cylinder: ",  baseName)
                sys.exit(6)            
        #End try
        forma=MAux.shape
        if (len(forma) >1):
            nLines=forma[0]
        else:
            nLines=1
        #End if len(forma)
        if nLines>2:
            try:
                print("Reading existingTime, cylinder: ",  baseName)
                existingTime=np.genfromtxt(outputName)[:, [0, 1]]
            except:
                print("Something went wrong with existingTime reading, cylinder: ",  baseName)
                print("Backing up to file .bak1 and Deleting the last line, cylinder: ",  baseName)
                copyFile(outputName, backupName+"1")
                deleteLastLine(outputName)
                try:
                    print("Reading the existingTime again, cylinder: ",  baseName)
                    existingTime=np.genfromtxt(outputName)[:, [0, 1]]
                except:
                    print("Something went double wrong with existingTime reading, cylinder: ",  baseName)
                    sys.exit(6)
            #End try
        else:
            print("The file only has ",  nLines, " lines for cylinder: ",  baseName)
            print("Backing up to file .bak2 and Deleting the complete file, cylinder: ",  baseName)
            copyFile(outputName, backupName+"2")
            os.remove(outputName) #Deleting the file if we only have 2 lines o less
        #End if nLines
    #End if os.path.isfile
    existingLastEnd=existingTime[-1, 1]
    existingLastIni=existingTime[-1, 0]    
    arrayTimeLast=arrayTimes[-1]
    print("existingLastIni: ", existingLastIni, " existingLastEnd: ",  existingLastEnd,  " arrayTimeLast: ",  arrayTimeLast)

    #Starting the postprocessing
    if (existingLastEnd<arrayTimeLast):
        if os.path.isfile(outputName):
            print("ExistingTime < arrayTimeLast. Coninuint with , cylinder: ",  baseName)
            print("Backing up to file .bak and Deleting the last line, cylinder: ",  baseName)
            copyFile(outputName, backupName)
            deleteLastLine(outputName)
            print("Reading the existingTime after deleting the last line, cylinder: ",  baseName)
            existingTime=np.genfromtxt(outputName)[:, [0, 1]]
        #End if os.path.isfile
        
        #Defining where to continue the analysis
        existingLastEnd=existingTime[-1, 1]
        existingLastIni=existingTime[-1, 0]
        startIndexTimes=int(np.searchsorted(arrayTimes,existingLastEnd)) #Starting from the  beginning of the last line
        print("Indexing the start of process")
        print("existingLastIni: ", existingLastIni, " existingLastEnd: ",  existingLastEnd)        

        #Initial reading of first data and ordering of angular samples ------------------------------
        jTimeCS=listTimesClean[startIndexTimes]
        jTimeNS=arrayTimes[startIndexTimes]
        print("Starting postprocess on time: ",  jTimeNS)
        dataNameS=os.path.join(volumetricDir, jTimeCS, iName)
        dataArrayS=np.genfromtxt(dataNameS, delimiter=None) #, usecols=(0, 1)) 
        dataArrayS=dataArrayS[rowIni:rowEnd, :]
            
        #Reading the coordinates and values, and converting to local polar coordinates
        xLocalS=dataArrayS[:, 0]-centresClean[ii, 0]
        yLocalS=dataArrayS[:, 1]-centresClean[ii, 1]
        zLocalS=dataArrayS[:, 2]
        thetaLocalS=np.arctan2(yLocalS, xLocalS)
        print("The size of thetaLocalS is: ", thetaLocalS.shape)
        rpLocalS=(np.sqrt(np.power(xLocalS, 2)+np.power(yLocalS, 2))-rCollector)/rCollector        
        UxLocalS=dataArrayS[:, 3]
        UyLocalS=dataArrayS[:, 4]
        UzLocalS=dataArrayS[:, 5]
        UrLocalS=np.multiply(UxLocalS, np.cos(thetaLocalS))+np.multiply(UyLocalS, np.sin(thetaLocalS))
        UthetaLocalS=-np.multiply(UxLocalS, np.sin(thetaLocalS))+np.multiply(UyLocalS, np.cos(thetaLocalS))

        #Checking the radial positions
        minRp=np.min(rpLocalS)
        maxRp=np.max(rpLocalS)
        diffRp=np.abs(maxRp-minRp)
        if diffRp>1E-6:                    
            print("Radial positions have a problem at S stage")
            print(rpLocalS)            
            sys.exit(7)
        #End if diffRp
                  
        #Sorting the arrays with the thetaLocal and checking the size of the read array
        zipAS=np.stack((thetaLocalS, rpLocalS, xLocalS, yLocalS, zLocalS, UrLocalS, UthetaLocalS), 1)
        sortingOrderS=np.argsort(zipAS[:, 0])
        zipSortS=zipAS[sortingOrderS]
        deltaThetaS=zipSortS[1:, 0]-zipSortS[0:-1, 0]
        deltaThetaS=np.reshape(np.insert(deltaThetaS, 0, np.amax(deltaThetaS)), (len(deltaThetaS)+1, 1))
        zipSortPlusS=np.append(zipSortS, deltaThetaS,axis=1)
        if (len(deltaThetaS)>NthetaDiv): #If OpenFOAM generated an extra sample point in its circles
            lineOutS=np.argmin(zipSortPlusS, axis=0)[-1]
            zipSortMinusS=np.delete(zipSortPlusS, lineOutS, axis=0)
        else:
            zipSortMinusS=zipSortPlusS
        #End if len>NthetaDiv
          
        #Integrating the fluxes towards the cylinder
        UrLocalS=zipSortMinusS[:, -3] #Third column from end to beginning is Ur
        rpLocalS=zipSortMinusS[:,1] #Second column (starting from 0 index) is rpLocal
        thetaLocalS=zipSortMinusS[:,0] #First column (starting from 0 index) is thetaLocal
        print("The size of thetaLocalS is: ", thetaLocalS.shape)
        NTheta=np.size(thetaLocalS)
        UrMClippedS=np.negative(np.clip(UrLocalS, -np.finfo('d').max, 0))
        kAreaS=deltaThetaS[0]*(rpLocalS[0]+1)*rCollector*deltaL[2]
        ##radialVolumeAccS=UrMClippedS*kAreaS*deltaT    
        
        #Writing the angles headers if it is the first line to write
        if (existingLastEnd == -1000):
            outputLine=np.reshape(np.insert(thetaLocalS,0,  [-1000, -1000]), (1, NTheta+2))
            print("The shape of outputLine is: ",  outputLine.shape)
            fOut=open(outputName, mode='ab')
            np.savetxt(fOut, outputLine)
            fOut.close()
        #End if existingLastEnd==-1000 (For printing the header
              

        #Finding the beginning and end of accumulation for the first time
        indexIni=startIndexTimes
        timeIni=arrayTimes[indexIni]
        timeEnd=np.floor(timeIni+1) #Tentative
        indexEnd=int(np.searchsorted(arrayTimes,timeEnd)) 
        indexEnd=np.minimum(indexEnd, Ntimes-1)
        timeEnd=arrayTimes[indexEnd]
        
        #Declaring the array for accumulation starting from zero
        radialVolumeAcc=np.zeros(UrMClippedS.shape)
        
        #If there is a single sampled time, i will consider it as a sample from Mean fields
        if (Ntimes == 1):
            deltaT=1
            radialVolumeAcc=radialVolumeAcc+(UrMClippedS+UrMClippedS)/2.0*kAreaS*deltaT
            
            #Writing the new volumetric integral into this time
            print("Writing a single data in the range: ", timeIni, " ",timeIni)
            outputLine=np.reshape(np.insert(radialVolumeAcc,0,  [timeIni, timeIni]), (1, NTheta+2))
            fOut=open(outputName, mode='ab')
            np.savetxt(fOut, outputLine)
            fOut.close()
            
        else:        
            
            for jj in range(startIndexTimes+1,Ntimes):           
                jTimeC=listTimesClean[jj]
                jTimeN=arrayTimes[jj]
                deltaT=jTimeN-jTimeNS
                
                
                dataName=os.path.join(volumetricDir, jTimeC, iName)
                #print(dataName)
                dataArray=np.genfromtxt(dataName, delimiter=None) #, usecols=(0, 1)) 
                dataArray=dataArray[rowIni:rowEnd, :]
                            
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

                #Checking the radial positions
                minRp=np.min(rpLocal)
                maxRp=np.max(rpLocal)
                diffRp=np.abs(maxRp-minRp)
                if diffRp>1E-6:                    
                    print("Radial positions have a problem")
                    print(rpLocal)
                    sys.exit(7)
                #End if diffRp
                          
                #Sorting the arrays with the thetaLocal and checking the size of the read array
                zipA=np.stack((thetaLocal, rpLocal, xLocal, yLocal, zLocal, UrLocal, UthetaLocal), 1)            
                sortingOrder=np.argsort(zipA[:, 0])
                zipSort=zipA[sortingOrder]
                deltaTheta=zipSort[1:, 0]-zipSort[0:-1, 0]
                deltaTheta=np.reshape(np.insert(deltaTheta, 0, np.amax(deltaTheta)), (len(deltaTheta)+1, 1))
                zipSortPlus=np.append(zipSort, deltaTheta,axis=1)        
                if (len(deltaTheta)>NthetaDiv): #If OpenFOAM generated an extra sample point in its circles
                    lineOut=np.argmin(zipSortPlus, axis=0)[-1]
                    zipSortMinus=np.delete(zipSortPlus, lineOut, axis=0)
                else:
                    zipSortMinus=zipSortPlus
                #End if len>NthetaDiv
                            
                #Integrating the fluxes towards the cylinder
                UrLocal=zipSortMinus[:, -3] #Third column from end to beginning is Ur
                rpLocal=zipSortMinus[:,1] #Second column (starting from 0 index) is rpLocal
                UrMClipped=np.negative(np.clip(UrLocal, -np.finfo('d').max, 0))
                kArea=deltaTheta[0]*(rpLocal[0]+1)*rCollector*deltaL[2]
                radialVolumeAcc=radialVolumeAcc+(UrMClipped+UrMClippedS)/2.0*kArea*deltaT

                #Final operations in the cycle
                
                if ((jj == indexEnd) or (jj == Ntimes-1)):
                    #Writing the new volumetric integral into this time
                    print("Writing data in the range: ", timeIni, " ",jTimeN)
                    outputLine=np.reshape(np.insert(radialVolumeAcc,0,  [timeIni, jTimeN]), (1, NTheta+2))
                    fOut=open(outputName, mode='ab')
                    np.savetxt(fOut, outputLine)
                    fOut.close()
                #End if jj==indexEnd
                
                if ((jj == indexEnd) and (jj != Ntimes-1)):
                    print("IndexEnd:", indexEnd)
                    print("Ntimes-1:", Ntimes-1)
                    print("Updating times from : ", timeIni, " ",timeEnd)
                    #Finding the beginning and end of accumulation for the first time
                    indexIni=indexEnd #Starting from the  beginning of the last line
                    timeIni=arrayTimes[indexIni]
                    timeEnd=np.floor(timeIni+1) #Tentative
                    indexEnd=int(np.searchsorted(arrayTimes,timeEnd))
                    indexEnd=np.minimum(indexEnd, Ntimes-1)
                    timeEnd=arrayTimes[indexEnd]
                    print("To                  : ", timeIni, " ",timeEnd)                
                    #Declaring the array for accumulation starting from zero
                    radialVolumeAcc=np.zeros(UrMClippedS.shape)
                #End if jj==indexEnd
                
                #Cyclic substitution
                UrMClippedS=UrMClipped
                jTimeNS=jTimeN
            #End for jj in range
        #End if (Ntimes == 1):
    #End if  existingLastEnd<arrayTimeLast  
    prueba=np.genfromtxt(outputName)
    print("Aqui tienes la prueba")
    print(prueba) 
#End of the for ii in range(0,NCylinders) .OR. def calcOneAngular function

#Declaring the number of parallel workers
num_of_workers = multiprocessing.cpu_count()
pool = multiprocessing.Pool(num_of_workers)

#Operating the calcOneAngular function in parallel
for data in pool.map(calcOneAngular, inputs):
    ii=data
        
#End of script
print("Script done")
    

