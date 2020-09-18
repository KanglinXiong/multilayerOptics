#########################################################
# This file is part of the MultilayerOptics package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

import os.path, math
import Optics, Misc

''' This module contains the API functions. '''

# The OptStructure object that will do the analysis.
optStruct = None

# The script name will be used in naming of result files.
scriptName = None

def buildStruct(scriptFilename):
    ''' build the optical structure, record the script name '''
    global optStruct, scriptName
    optStruct = Optics.OptStructure(scriptFilename)
    scriptName = os.path.splitext(scriptFilename)[0]

def setMaterialAbsorption(material, a):
    ''' The unit for absorption coefficient is 1/cm. '''
    optStruct.setMaterialGain(material, 0 - a)

def searchLaserMode(gainMedium, g0, w0, g1, g2, w1, w2):
    return(optStruct.searchMode(gainMedium, g0, w0, g1, g2, w1, w2))

def saveReflect(w1, w2, dw):
    ''' [[x, r], ...] '''
    filename = scriptName + " "  + "reflectivity.txt"
    reflectFile = open(Misc.getResultfilename(filename), 'w')
    for r in optStruct.getReflectSpectrum(w1, w2, dw):
        reflectFile.write(str(r[0])+"\t"+str(r[1])+"\n")
    reflectFile.close()

def saveField(w, dx):
    ''' [[x, |E|**2], ...] '''
    filename = scriptName + " "  + "field.txt"
    fieldFile = open(Misc.getResultfilename(filename), 'w')
    for f in optStruct.getFieldDistribution(w, dx):
        fieldFile.write(str(f[0])+"\t"+str(f[1])+"\n")
    fieldFile.close()

def saveIndex():
    ''' [[x, n], ...] '''
    filename = scriptName + " "  + "index.txt"
    indexFile = open(Misc.getResultfilename(filename), 'w')
    for n in optStruct.getRefractiveIndexDistribution():
        indexFile.write(str(n[0])+"\t"+str(n[1])+"\n")
    indexFile.close()


def run(scriptFilename):
    ''' Parse the script file and evaluate the commands behind the #> prompt. '''
    # Build the OptStructure using layers defined in the script.
    buildStruct(scriptFilename)
    print("Start to execute the script ", scriptName)
    # Execute the commands in the script.
    fileObj = open(Misc.getScriptfilename(scriptFilename), mode="rt", newline='\n')
    for line in fileObj:
            line = line.strip()
            if(len(line) < 3 or not line.startswith("#>")):
                continue
            print(line)
            eval(line.strip("#>"))
    print("Finished executing the script ", scriptName)
    print("Press <Enter> to exit.")
    input()


def computeOpticalLength(w, matNameAndThickList):
    '''
    A quick method to calculate optical length of a stack of layers at  given wavelength.
    The 2nd parameter is [mat1, d1, mat2, d2, ..., matN, dN].
    The angle is not considered yet.
    '''
    count = round(len(matNameAndThickList)/2)
    if(count*2 < len(matNameAndThickList)*0.99 or count*2 > len(matNameAndThickList)*1.01):
        print(matNameAndThickList, "is invalid input.")
        exit()
    length = 0.0
    for i in range(count):
        [matName, thickness] = [matNameAndThickList[2*i], matNameAndThickList[2*i+1]]
        length += (Optics.OptMaterial(matName).getNk(w).real)*thickness
    return(length)


if __name__ == "__main__":
    run(scriptFilename = "2016-05-04 blue VCSEL.txt")



