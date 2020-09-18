#########################################################
# This file is part of the MultilayerOptics package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################


import os, os.path, re

''' This module manage the locations of files, and expand graded layers. '''

def getNkfilename(filename):
    ''' The nk data is in a sub directory. '''
    directory = os.path.dirname(os.path.realpath(__file__))
    return(os.path.join(directory, "nk", filename))

def getScriptfilename(filename):
    ''' The script is in a sub directory. '''
    directory = os.path.dirname(os.path.realpath(__file__))
    return(os.path.join(directory, "script", filename))

def getResultfilename(filename):
    ''' The results will be in a sub directory. '''
    dir = os.path.dirname(os.path.realpath(__file__))
    return(os.path.join(dir, "result", filename))

def interpolate2NkFiles(x0, x1, x2, filename1, filename2):
    ''' Given nk files of Ax1BC and Ax2BC, get nk data of Ax0BC. '''
    # read nk files of Ax1BC and Ax2BC
    def readNk(filename):
        ''' the nk data is [[w, n, k], ...] with increasing wavelength. '''
        nk = []
        fileObj = open(getNkfilename(filename), mode="rt", newline='\n')
        fileObj.readline()
        for line in fileObj:
            tmp = line.replace("\t", " ").replace(",", " ").split(" ")
            nk.append([float(tmp[0]), complex(float(tmp[1]), float(tmp[2]))])
        return(nk)
    
    # interpolate wavelength of Ax1BC and Ax2BC
    def getNkByW(w, nk):
        ''' nk is [[w, complex(n, k)], ...]. '''
        if(w <= nk[0][0]):
            return(nk[0][1])
        elif(w >= nk[-1][0]):
            return(nk[-1][1])
        i = 0
        while(i < len(nk) and w > nk[i][0]):
            i += 1
        if(i == len(nk) - 1):
            return(nk[-1][1])
        [d1, d2] = [w - nk[i][0], nk[i+1][0]- w]
        return((nk[i][1]*d2 + nk[i+1][1]*d1)/(d1+d2))
    
    # interpolate alloy composition for Ax0BC
    [nk0, nk1, nk2] = [[], readNk(filename1), readNk(filename2)]
    wList = []
    for i in range(len(nk1)):
        wList.append(nk1[i][0])
    for i in range(len(nk2)):
        if not (nk2[i][0] in wList):
            wList.append(nk2[i][0])
    wList.sort(reverse = False)   
    for w in wList:
            y = ((x0 - x1)*getNkByW(w, nk2) + (x2 - x0)*getNkByW(w, nk1))/(x2-x1)
            nk0.append([w, y])

    # the nk data for alloy Ax0BC
    return(nk0)

def getTernaryAlloyNk(alloyName = "Al25%GaAs"):
    '''
    Find nk files of alloys and interplotate the nk.
    '''
    directory = os.path.dirname(os.path.realpath(__file__))
    directory = os.path.join(directory, "nk")
    filenameList = os.listdir(directory)

    digits = re.findall("[0-9]{2}", alloyName)
    if len(digits) != 1:
        print("Error", __file__)
    digits = digits[0]

    pattern = alloyName.replace(digits, "[0-9]{2}")
    regex = re.compile(pattern + ".txt")
    [matchList, digitsList] = [[], []]
    for filename in filenameList:
        rlt = regex.match(filename)
        if(rlt != None):
            matchList.append(rlt.string)
            digitsList.append(re.findall("[0-9]{2}", rlt.string)[0])
    if(len(matchList) < 2):
        print("Error", __file__)

    digitsList = list(map(int, digitsList))
    digitsList.sort()
    digits = int(digits)
    for i in range(len(digitsList)-1):
        if (digitsList[i] - digits)*(digitsList[i+1] - digits) < 0:
            break
    if (digitsList[i] - digits)*(digitsList[i+1] - digits) >= 0:
        print("Error", __file__)

    [filename1, filename2] = [matchList[i], matchList[i+1]]
    [x0, x1, x2] = [digits, digitsList[i], digitsList[i+1]]
    return(interpolate2NkFiles(x0, x1, x2, filename1, filename2))


def discretizeGradedLayer(cmd = "[Al60%02%GaAs 120.0 21]"):
    '''
    the input command is [Ax1%x2BC% thickness numberOflayers].
    The output is a line of string for OptStructure.setStack.
    '''
    eleList = cmd.strip("[]").split()
    if len(eleList) != 3:
        print("Error", __file__)
    [name, d, n] = [eleList[0], float(eleList[1]), int(eleList[2])]
    if n <= 2:
        print("Error", __file__)
    [x1, x2] = list(map(int, re.findall("[0-9]{2}", name)))
    dx = (x2 - x1)/(n - 1.0)
    xStrList = list(map(lambda i: str(round(x1 + dx*i)).zfill(2) + "%", range(n)))
    labelList = re.split("[0-9]{2}%[0-9]{2}%", name)
    strList = list(map(lambda i: labelList[0] + xStrList[i] + labelList[1], range(n)))
    string = ""
    for i in range(n - 1):
        string = string + strList[i] + " " + str(d/float(n)) + ","
    string = string + strList[-1] + " " + str(d/float(n))
    print("graded layer ", cmd, "will be expanded as \n", string)
    return(string)


if __name__ == "__main__":
    print(getTernaryAlloyNk("Al19%GaAs"))

