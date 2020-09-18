#########################################################
# This file is part of the MultilayerOptics package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################


# This module is to analyze optical properties of a stack of layers.
# The numerical accuracy is not limited by cmath.exp() or cmath.abs().
# The data [[wavelength, n, k], ...] is linearly interpolated.
# The file name of nk should only have [a-Z0-9%], without other symbols.

import math, cmath, copy, os.path, re
import Solver, Misc

'''
Units:         length in nm, gain in 1/cm
Axis:          0 ----> x, that is, from left to right
Structure:     layer 1 | layer 2 | ... | layer n
Waves:         <-- left traveling, right traveling --> 
'''

class OptMaterial:
    '''
    The material used by layer.
    '''
    def __init__(self, name = None):
        self.name = name
        self.nk = []
        self.gain = None
        if(name != None):
            self.setNk(name + ".txt")
        
    def getNk(self, wavelength):
        '''
        The {wavelength, n, k} is linearly interpolated.
        '''
        nk = complex(1, 0)
        if (wavelength <= 0.0):
            return nk
        if(wavelength >= self.nk[-1][0]):
            nk = self.nk[-1][1]
        elif(wavelength <= self.nk[0][0]):
            nk = self.nk[0][1]
        else:
            i = 0
            while(i < len(self.nk) and wavelength > self.nk[i][0]):
                i += 1
            if(i == len(self.nk) - 1):
                nk = self.nk[-1][1]
            else:
                d1 = wavelength - self.nk[i][0]
                d2 = self.nk[i+1][0]- wavelength
                nk = (self.nk[i][1]*d2 + self.nk[i+1][1]*d1)/(d1+d2)
        if(self.gain != None):
            imag = (0.0 - self.gain) * (wavelength * 1e-7) / (4.0 * math.pi)
            return complex(nk.real, imag)
        else:
            return nk

    def setNk(self, filename):
        '''
        The index file is in the format of {wavelength, n, k} with increasing wavelength. 
        '''
        # need to interplotate nk of alloys
        if not os.path.exists(Misc.getNkfilename(filename)):
            self.nk = Misc.getTernaryAlloyNk(filename.replace(".txt", "").strip())
            return
        # the nk files exists
        self.nk = []
        fileObj = open(Misc.getNkfilename(filename), mode="rt", newline='\n')
        fileObj.readline()
        for line in fileObj:
            tmp = line.replace("\t", " ").replace(",", " ").split(" ")
            self.nk.append([float(tmp[0]), complex(float(tmp[1]), float(tmp[2]))])

    def setGain(self, g):
        '''
        set gain to g (1/cm)
        '''
        self.gain = g

    def setAbsorption(self, a):
        '''
        set gain to -a (1/cm)
        '''
        self.gain = 0.0 - a


class OptLayer:
    '''
    The amplitude and phase values of electrical waves are defined at the center of the layer.
    '''
    def __init__(self, material = None):
        #OptMaterial instance
        self.material = material 
        self.thickness = 0.0
        self.center = 0.0
        self.angle = 0.0
        self.wavelength = 520.0
        #[complex amplitude, phase] of left-traveling wave
        self.El = [complex(1,0), 0]
        #[complex amplitude, phase] of right-traveling wave
        self.Er = [complex(1,0), 0]

    def setLayer(self, material, thickness):
        self.material = material
        self.thickness = thickness

    def setWavelegnth(self, wavelength):
        self.wavelength = wavelength

    def setAngle(self, angle):
        self.angle = angle

    def getThickness(self):
        return self.thickness

    def getNk(self):
        return self.material.getNk(self.wavelength)

    def getAngle(self):
        return self.angle

    def getElec(self, displace = 0):
        '''
        n affects phase, k affects amplitude.
        Displacement is -thickness/2 for left, 0 for center, thickness/2 for right.
        '''
        nk = self.material.getNk(self.wavelength)
        waveVec = 2 * math.pi * nk / self.wavelength
        complexPhase = waveVec * (displace / math.cos(self.angle)) * complex(0, 1)
        El = copy.deepcopy(self.El)
        Er = copy.deepcopy(self.Er)
        El[0] = El[0] / cmath.exp(complexPhase)
        Er[0] = Er[0] * cmath.exp(complexPhase)
        El[1] = El[1] - complexPhase.imag
        Er[1] = Er[1] + complexPhase.imag
        return(El, Er)

    def setElec(self, El, Er, displace = 0):
        '''
        Set El and Er at the layer center by values at given displacement.
        '''
        nk = self.material.getNk(self.wavelength)
        waveVec = 2 * math.pi * nk / self.wavelength
        complexPhase = waveVec * ((0.0 - displace) / math.cos(self.angle)) * complex(0, 1)
        self.El[0] = El[0] / cmath.exp(complexPhase)
        self.Er[0] = Er[0] * cmath.exp(complexPhase)
        self.El[1] = El[1] - complexPhase.imag
        self.Er[1] = Er[1] + complexPhase.imag

    def getRightSideElec(self):
        return self.getElec(self.thickness / 2.0)

    def getLeftSideElec(self):
        return self.getElec(0 - self.thickness / 2.0)


class OptStructure():
    '''
    From left to Right, the stack is layer 1 | layer 2 | ... |layer n.
    The script is {mat1 d1, mat2 d2, ... } * n, or {mat, d} from left to right.
    '''
    def __init__(self, filename = ''):
        self.wavelength = 520
        self.S = True
        self.P = False
        self.polarization = self.S
        # [[mat1, d1], [mat2, d2], ...]
        self.stack = []
        # [optLayer1, optLayer2, ...]
        self.struct = []
        self.matDict = None
        if(filename):
            self.setStack(filename)
            self.setStruct()

    def setWavelength(self, wavelength = 520):
        '''
        Set wavelength of each OptLayer instance.
        '''
        self.wavelength = wavelength
        for i in range(len(self.struct)):
            self.struct[i].wavelength = self.wavelength

    def setStruct(self):
        '''
        The self.struct is a list of OptLayer instances, [optLayer1, optLayer2, ...].
        '''
        tmpMatList = []
        for lay in self.stack:
            if(not lay[0] in tmpMatList):
                tmpMatList.append(lay[0])
        tmpMatDict = []
        for mat in tmpMatList:
            tmpMatDict.append([mat, OptMaterial(mat)])
        self.matDict = dict(tmpMatDict)
        self.struct = []
        for lay in self.stack:
            self.struct.append(OptLayer(self.matDict[lay[0]]))
            self.struct[-1].thickness = lay[1]

    def setStack(self, filename):
        '''
        Parse the structure file to set the self.stack.
        The self.stack is a list of material and thickness, [[mat1, d1], [mat2, d2], ...].
        '''
        fileObj = open(Misc.getScriptfilename(filename), mode="rt", newline='\n')
        for line in fileObj:
            line = line.strip()
            if((not len(line)) or ("#" in line)):
                continue
            line = line.expandtabs(1)
            while(" " * 2 in line):
                line = line.replace(" " * 2, " ")
            layers = line.partition("{")[2].partition("}")[0]
            repeat = line.partition("*")[2].strip()
            # discretize graded layers like [Al60%02%GaAs 120.0 21], if any.
            gradedLayerList = re.findall("\[.*?\]", layers)
            for gradedLayer in gradedLayerList:
                layers = layers.replace(gradedLayer, Misc.discretizeGradedLayer(gradedLayer))
            # parse all the layers, get [[mat, d], ...]
            tmpStack = []
            if(len(repeat)):
                repeat = int(repeat)
            else:
                repeat = 1
            for lay in layers.split(","):
                lay = lay.strip().split(" ")
                lay[1] = float(lay[1])
                tmpStack.append(lay)
            for lay in (tmpStack * repeat):
                self.stack.append(lay)

    def getOpticalLength(self, wavelength = 520):
        ''' sum of n*d of all OptLayer '''
        self.setWavelength(wavelength)
        optLen = 0.0
        for i in range(len(self.struct)):
            optLen += self.struct[i].getNk() * self.struct[i].getThickness()
        return(optLen.real)

    def getInterfaceMatrix(self, nk1, theta1, nk2, theta2):
        '''
        Get transfer matrix for interface, by the boundary requirements of Maxwell's equations.
        Label 1 for input side, 2 for output side.
        The angle theta1 and theta2 should be calculated somewhere else by Snell's law.
        '''
        if(self.polarization == self.S):
            m1122 = cmath.cos(theta2) * nk2 - cmath.cos(theta1) * nk1
            m1221 = cmath.cos(theta2) * nk2 + cmath.cos(theta1) * nk1
            m0 = 2.0 * cmath.cos(theta2) * nk2
            return [[m1122/m0, m1221/m0], [m1221/m0, m1122/m0]]
        else:
            m1122 = cmath.cos(theta2) * nk1 - cmath.cos(theta1) * nk2
            m1221 = cmath.cos(theta2) * nk1 + cmath.cos(theta1) * nk2
            m0 = 2.0 * cmath.cos(theta2) * nk2
            return [[m1122/m0, m1221/m0], [m1221/m0, m1122/m0]] 

    def calcElecOfRightLayByLeftLay(self, leftLay, rightLay):
        '''
        The phase is calculated in the accumulated way.
        ''' 
        m = self.getInterfaceMatrix(leftLay.getNk(), leftLay.getAngle(),
                                   rightLay.getNk(), rightLay.getAngle())
        [inEl, inEr] = leftLay.getRightSideElec()
        # amplitudes of waves
        [inWl, inWr] = [inEl[0], inEr[0]]
        [outWl, outWr] = [m[0][0] * inWr + m[0][1] * inWl,
                          m[1][0] * inWr + m[1][1] * inWl]
        [outEl, outEr] = [[outWl, cmath.phase(outWl)],
                          [outWr, cmath.phase(outWr)]]
        while(outEl[1] > inEl[1] and inEl[0] != 0.0):
            outEl[1] -= 2.0 * math.pi
        while(outEr[1] < inEr[1] and inEr[0] != 0.0):
            outEr[1] += 2.0 * math.pi
        rightLay.setElec(outEl, outEr, 0.0 - rightLay.thickness/2.0)

    def calcElecOfLeftLayByRightLay(self, rightLay, leftLay):
        m = self.getInterfaceMatrix(rightLay.getNk(), rightLay.getAngle(),
                                 leftLay.getNk(), leftLay.getAngle())
        [inEl, inEr] = rightLay.getLeftSideElec()
        # amplitudes of waves
        [inWl, inWr] = [inEl[0], inEr[0]]                             
        [outWl, outWr] = [m[1][0] * inWl + m[1][1] * inWr,
                          m[0][0] * inWl + m[0][1] * inWr]
        [outEl, outEr] = [[outWl, cmath.phase(outWl)],
                          [outWr, cmath.phase(outWr)]]
        while(outEl[1] < inEl[1] and inEl[0] != 0.0):
            outEl[1] += 2.0 * math.pi
        while(outEr[1] > inEr[1] and inEr[0] != 0.0):
            outEr[1] -= 2.0 * math.pi
        leftLay.setElec(outEl, outEr, leftLay.thickness/2.0)

    def getLeftIncidentByRightExit(self, wavelength = 520):
        '''
        Given the right most, compute the left most.
        '''
        self.setWavelength(wavelength)
        #set right most
        self.struct[-1].setElec([0, 0], [1e0, 0])
        i = -1
        while(i > 0 - len(self.struct)):
            self.calcElecOfLeftLayByRightLay(self.struct[i], self.struct[i - 1])
            i -= 1
        #get left most
        [ref, inc] = self.struct[0].getRightSideElec()
        return [inc, ref]

    def getRightIncidentByLeftExit(self, wavelength = 520):
        '''
        Given the left most, compute the right most.
        '''
        self.setWavelength(wavelength)
        #set left most
        self.struct[0].setElec([1e0, 0], [0, 0])
        i = 0
        while(i < len(self.struct) - 1):
            self.calcElecOfRightLayByLeftLay(self.struct[i], self.struct[i + 1])
            i += 1
        #get right most
        [inc, ref] = self.struct[-1].getLeftSideElec()
        return [inc, ref]

    def getReflectOfLeftSurface(self, wavelength = 520):
        '''
        Get the reflectivity of the left surface of the structure.
        '''
        [inc, ref] = self.getLeftIncidentByRightExit(wavelength)
        if(abs(inc[0]) != 0.0):
           return (abs(ref[0]/inc[0])**2.0, ref[1] - inc[1])
        else:
           return ("Error!")

    def getReflectOfRightSurface(self, wavelength = 520):
        '''
        Get the reflectivity of the right surface of the structure.
        '''
        [inc, ref] = self.getRightIncidentByLeftExit(wavelength)
        if(abs(inc[0]) != 0.0):
           return (abs(ref[0]/inc[0])**2.0, ref[1] - inc[1])
        else:
           return ("Error!")

    def getReflectSpectrum(self, w1 = 500, w2 = 600, dw = 1, surface = "R"):
        '''
        Get reflectivity over a given wavelength range.
        '''
        reflectivity = []
        w = w1
        if(surface == "L"):
            while(w <= w2):
                reflectivity.append([w, self.getReflectOfLeftSurface(w)[0]])
                w = w + dw
        else:
            while(w <= w2):
                reflectivity.append([w, self.getReflectOfRightSurface(w)[0]])
                w = w + dw
        return reflectivity

    def getFieldIntensity(self, w = 520, loc = 0):
        '''
        Get the electrical field intensity within the OptStructure by a location. 
        If to get the energy intensity, use epsr ~ n**2, I ~ eps0*epsr*Elec**2.
        '''
        #The member function self.getReflectSpectrum should have been called.
        x = 0.0
        layer = None
        i = 0
        while(i < len(self.struct)):
            layer = self.struct[i]
            if(loc <= 0.0):
                x = 0.0 - layer.thickness / 2.0
                break
            elif(x <= loc and x + layer.thickness > loc):
                x = loc - (x + layer.thickness / 2.0)
                break
            x += layer.thickness
            i += 1
        if(i == len(self.struct) and x < loc):
            x = layer.thickness / 2.0
        [El, Er] = layer.getElec(x)
        return (abs(El[0] + Er[0])**2)

    def getFieldDistribution(self, w = 520, x1 = -1e5, x2 = 1e5, dx = 1, surface = "R"):
        '''
        Build the field and return the distribution [[x, intensity], ...].
        '''
        # Build the field in the OptStructure by call the number function.
        self.getReflectSpectrum(w, w, 1e1, surface)
        xMin = max(0.0, x1)
        xMax = 0.0
        for layer in self.struct:
            xMax += layer.thickness
        xMax = min(xMax, x2)
        field = []
        x = xMin
        while(x <= xMax):
            field.append([x, self.getFieldIntensity(w, x)])
            x += dx
        return field

    def getRefractiveIndexDistribution(self):
        ''' Return the [[x, n], ] distribution. '''
        nkList = []
        x = 0.0
        for layer in self.struct:
            if(layer.thickness == 0.0):
                continue
            dx = min(1e-3, layer.thickness * 1e-3)
            nkList.append([x + dx, layer.getNk().real])
            dx = max(layer.thickness - 1e-3, layer.thickness * (1.0 - 1e-3))
            nkList.append([x + dx, layer.getNk().real])
            x += layer.thickness
        return nkList

    def setMaterialGain(self, mat = "", gain = None):
        if(not mat in self.matDict):
            return None
        self.matDict[mat].setGain(gain)

    def searchMode(self, gainMedium = "QW", g0 = 1500, w0 = 520,\
                   g1 = 50, g2 = 5000, w1 = 100, w2 = 1600):
        '''
        Find mode and threshold gain by tuning material gain.
        When no incident light is needed, the cavity emits laser.
        refer to Analysis of multielement semiconductor lasers,
        Journal of Applied Physics 54, 2962 (1983)
        '''
        def getIncidentWave(g, w):
            ''' local function for Newton's method '''
            self.setMaterialGain(gainMedium, g)
            return abs(self.getRightIncidentByLeftExit(w)[0][0])
        # Find root for function getIncidentWave(g, w) within given tolerance.
        nt = Solver.Newton(getIncidentWave, g0, w0, g1, g2, w1, w2)
        print("begin mode searching...")
        print("iter, residual, gain(1/cm), wavelength(nm)")
        return (nt.run())


