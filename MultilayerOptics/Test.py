#########################################################
# This file is part of the MultilayerOptics package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################


import App

''' an example '''

w0 = 940

App.run(scriptFilename = "2020-07-18 tj vcsel.txt")

print("optical length: ", App.optStruct.getOpticalLength(w0))

print("refractive index:")
n1 = Optics.OptMaterial("Al15%GaAs").getNk(w0).real
n2 = Optics.OptMaterial("Al90%GaAs").getNk(w0).real
n3 = Optics.OptMaterial("Al20%GaAs").getNk(w0).real
print("Al15%GaAs", "Al90%GaAs", "Al20%GaAs")
print(n1, n2, n3)
