#########################################################
# This file is part of the MultilayerOptics package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

import random

''' This module is to solve f(x, y) == 0 by Newton's method. '''

class Newton:
    '''
    Newton's method to find root for f(x, y) == 0
    x(n+1) = x(n) - f(x(n))/f'(x(n)
    By default, f(x, y) is incidentLight(gain, wavelength)
    '''
    def __init__(self, f, x0, y0, x1, x2, y1, y2):
        '''
        x0 and y0 are initial values. [x1, x2] and [y1, y2] are ranges.
        '''
        self.f = f
        self.x = x0
        self.y = y0
        self.xmin = x1
        self.xmax = x2
        self.ymin = y1
        self.ymax = y2
        self.tolerance = 1e-3

    def dfdx(self):
        return (self.f(self.x + 1e-6, self.y)\
                - self.f(self.x, self.y))/1e-6

    def dfdy(self):
        return (self.f(self.x, self.y + 1e-6)\
                - self.f(self.x, self.y))/1e-6

    def run(self, maxIteration = 1e2):
        ''' Search with given tolerance and max iterations. '''
        count = 0
        while(count < maxIteration):
            count += 1
            # update x
            if(self.dfdx() != 0 ):
                dx = - self.f(self.x, self.y)/self.dfdx()
                if(self.x + dx < self.xmin or self.x + dx > self.xmax):
                    self.x = random.uniform(self.xmin, self.xmax)
                else:
                    self.x += dx
            else:
                self.x = random.uniform(self.xmin, self.xmax)
            # root is found
            if(abs(self.f(self.x, self.y)) < self.tolerance):
                break
            
            # update y
            if(self.dfdy() != 0): 
                dy = - self.f(self.x, self.y)/self.dfdy()
                if(self.y + dy < self.ymin):
                    self.y = self.ymin
                elif(self.y + dy > self.ymax):
                    self.y = self.ymax
                else:
                    self.y += dy
            else:
                self.y = random.uniform(self.xmin, self.xmax)
            # root is found
            if(abs(self.f(self.x, self.y)) < self.tolerance):
                break

            # print process
            print(count, self.f(self.x, self.y), self.x, self.y)
        return (self.f(self.x, self.y), self.x, self.y)
            

    
        
        
        
