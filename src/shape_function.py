from math import sin, pi, log10
import numpy as np
import sys

class ShapeFunction:
    '''Hicks-Henne or NACA shape function'''

    def __init__(self, x, domainidx=None):
        self.x = x
        npt = x.shape[0]
        self.shape = np.zeros((npt))

        # Optimization limits
        if domainidx == None:
            self.domainidx = [0, npt-1]
        elif domainidx[1] - domainidx[0] < 0:
            sys.stderr.write("Shape function domain must be positive.\n")
            sys.exit(1)
        elif (domainidx[0] < 0) or (domainidx[1] >= npt):
            sys.stderr.write("Shape function domain is not valid.\n")
            sys.exit(1)
        self.lidx = domainidx[0]
        self.ridx = domainidx[1]

        # Rescale x so that modified portion is in the range [0,1]
        # Since the shape function is supposed to be identically 0 at xscale=0 and 1,
        #   leave those points out of xscale.
        xl = self.x[self.lidx]
        xr = self.x[self.ridx]
        self.xscale = (self.x[self.lidx+1:self.ridx] - xl)/(xr - xl)

    def createHicksHenne(self, st, t1, t2):
        '''Creates Hicks-Henne shape function from parameters
        
        Inputs:
            st: strength / scaling factor. st = 1 results in bump with height 1
            t1: bump location
            t2: bump width

        Returns:
            None, but computes self.shape
        '''

        # Enforce side constraints on bumps
        if t1 <= 0.:
            t1 = 0.001
        if t1 >= 1.:
            t1 = 0.099
        if t2 <= 0.:
            t2 = 0.001

        # Create shape function
        power1 = log10(0.5)/log10(t1)
        self.shape[self.lidx+1:self.ridx] = st*np.sin(pi*self.xscale**power1)**t2

    def createNACA(self, shapenum):
        '''Creates NACA shape function given a shape index.'''

        # First mode
        if shapenum == 1:
            self.shape[self.lidx+1:self.ridx] = np.sqrt(self.xscale) - self.xscale
        # Whole powered modes
        elif shapenum%2 == 0:
            power = float(shapenum)/2.
            self.shape[self.lidx+1:self.ridx] = self.xscale**power(1. - self.xscale)
        # Fractional powered modes
        else:
            power1 = float(shapenum-1)/2. + 2.
            power2 = power1 - 1.
            self.shape[self.lidx+1:self.ridx] = self.xscale**power1 - self.xscale**power2

        # Normalize
        maxval = np.max(abs(self.shape[self.lidx+1:self.ridx]))
        self.shape /= maxval
