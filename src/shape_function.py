from math import pi, log10
import numpy as np

class ShapeFunction:
    '''Hicks-Henne shape function'''

    def __init__(self, x, domain=[0.,1.]):
        self.x = x
        npt = x.shape[0]
        self.shape = np.zeros((npt))
        self.domain = domain            # Optimization limits

    def createShape(self, st, t1, t2):
        '''Creates Hicks-Henne shape function from parameters
        
        Inputs:
            st: strength / scaling factor. st = 1 results in bump with height 1
            t1: bump location
            t2: bump width

        Returns:
            None, but computes self.shape
        '''

        # Enforce side constraints on bumps
        if t1 <= self.domain[0]:
            t1 = self.domain[0] + 0.001
        if t1 >= self.domain[1]:
            t1 = self.domain[1] - 0.001
        if t2 <= 0.:
            t2 = 0.001

        # Create shape function
        power1 = log10(0.5)/log10(t1)
        self.shape = st*np.sin(pi*self.x**power1)**t2
