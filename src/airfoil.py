import sys
import numpy as np
import xml.etree.ElementTree as ET
import libxfoil_wrap as lxf
from libxfoil import xfoil_data_group
from shape_function import ShapeFunction

class Airfoil:
    def __init__(self, npt=0):
        self.x = np.zeros((npt))
        self.y = np.zeros((npt))
        self.xdg = xfoil_data_group()
        lxf.xfoil_init(self.xdg)
        self.xle = None
        self.yle = None

    def numPoints(self):
        return self.x.shape[0]

    def smoothPaneling(self, xfoilpanelingsettings):
        '''Smooths airfoil using libxfoil routine
        
        Inputs:
            xfoilpanelingsettings: what the name implies
        
        Returns:
            retval: True on success, False on error
            errmsg: Error message text
        '''

        errmsg = "Success"
        if self.numPoints() < 3:
            errmsg = "Cannot smooth paneling. Airfoil must have at least 3 points."
            sys.stderr.write(errmsg+"\n")
            return False, errmsg

        lxf.xfoil_set_buffer_airfoil(self.xdg, self.x, self.y, self.numPoints())
        lxf.xfoil_set_paneling(self.xdg, xfoilpanelingsettings.toXfoilGeomOpts())
        if lxf.xfoil_smooth_paneling(self.xdg) != 0:
            errmsg = "Error in xfoil_smooth_paneling. Please report this issue."
            sys.stderr.write(errmsg+"\n")
            return False, errmsg

        x, y, stat = lxf.xfoil_get_current_airfoil(self.xdg, xfoilpanelingsettings.value("npan"))
        if stat != 0:
            errmsg = "Error in xfoil_get_current_airfoil. Please report this issue."
            sys.stderr.write(errmsg+"\n")
            return False, errmsg
        self.x = np.array(x)
        self.y = np.array(y)

        return True, errmsg

    def findLE(self):
        '''Finds leading edge
        
        Returns:
            retval: True on success, False on error
            errmsg: message describing the error
        '''
        retval = True
        errmsg = "Success"
        if self.numPoints() < 3:
            errmsg = "Airfoil must have at least 3 points."
            return False, errmsg
        s, xs, ys = lxf.xfoil_spline_coordinates(self.x, self.y, self.numPoints())
        _, self.xle, self.yle = lxf.xfoil_lefind(self.x, self.y, s, xs, ys, self.numPoints())

        return retval, errmsg

    def chord(self):
        '''Computes chord. findLE must be called first.

        Returns:
            chord: airfoil chord
            stat: True on success, False on error
            errmsg: message describing the error
        '''
        errmsg = "Success"
        if self.numPoints() < 2:
            errmsg = "Airfoil must have at least 2 points."
            return 0., False, errmsg
        if self.xle is None or self.yle is None:
            errmsg = "findLE must be called first."
            return 0., False, errmsg
        xte = 0.5*(self.x[0] + self.x[self.numPoints()-1])

        return xte - self.xle, True, errmsg

    def translate(self, dx, dy):
        '''Translates airfoil'''
        self.x += dx
        self.y += dy
        self.xle += dx
        self.yle += dy

    def scale(self, xscale, yscale):
        '''Scales airfoil by factors in x and y'''
        self.x *= xscale
        self.y *= yscale
        self.xle *= xscale
        self.yle *= yscale


class SeedAirfoil(Airfoil):

    npointside = 100

    def __init__(self):
        super(SeedAirfoil, self).__init__()
        self.leidx = None
        self.addpoint_loc = None
        self.xt = np.zeros((0))
        self.yt = np.zeros((0))
        self.xb = np.zeros((0))
        self.yb = np.zeros((0))
        self.source_data = {"source": None, "camber": None, "xcamber": None,
                            "thickness": None, "designation": None, "file": None}
        self.domaint = [0.,1.]
        self.domainb = [0.,1.]
        self.domainidxt = None
        self.domainidxb = None
        self.shapest = []
        self.shapesb = []

    def readFromFile(self, fname):
        """Reads airfoil from file.

        Inputs
            fname: file name or path
        Returns
            retval: True on success, False on error
            errmsg: message describing the error
        """
        
        errmsg = "Success"
        print("Reading airfoil file {:s}".format(fname))
        try:
            f = open(fname)
        except IOError:
            errmsg = "Unable to read {:s}.".format(fname)
            sys.stderr.write(errmsg+"\n")
            return False, errmsg

        # Determine if file has a label
        labeled = False
        splitline = f.readline().split()
        if len(splitline) < 2:
            labeled = True
        else:
            try:
                x = float(splitline[0])
                y = float(splitline[1])
            except ValueError:
                labeled = True
        if labeled:
            print("Labeled coordinates file.")
        else:
            print("Plain coordinates file.")
        f.seek(0)

        # Read coordinates
        linenum = 0
        x = []
        y = []
        for line in f:
            if not labeled or linenum > 0:
                splitline = line.split()
                if len(splitline) < 2:
                    f.close()
                    errmsg = "Not enough data in line."
                    sys.stderr.write(errmsg+"\n")
                    return False, errmsg
                try:
                    x.append(float(splitline[0]))
                    y.append(float(splitline[1]))
                except ValueError:
                    f.close()
                    errmsg = "Unable to convert to float."
                    sys.stderr.write(errmsg+"\n")
                    return False, errmsg
            linenum += 1
        if len(x) > 2:
            self.x = np.array(x)
            self.y = np.array(y)
        else:
            f.close()
            errmsg = "Less than 3 vertices in coordinates file."
            sys.stderr.write(errmsg+"\n")
            return False, errmsg
        f.close()

        # FIXME: make sure ordering is counterclockwise

        self.source_data["source"] = "file"
        self.source_data["file"] = fname
        return True, errmsg

    def generate4Digit(self, camber, xcamber, thickness):
        x, y, _ = lxf.naca_4_digit(camber, xcamber, thickness, Airfoil.npointside)
        self.x = np.array(x)
        self.y = np.array(y)
        self.source_data["source"] = "4-digit"
        self.source_data["camber"] = camber
        self.source_data["xcamber"] = xcamber
        self.source_data["thickness"] = thickness

    def generate5Digit(self, designation):
        x, y, _, stat = lxf.naca_5_digit(designation, Airfoil.npointside)
        if stat != 0:
            return False
        else:
            self.x = np.array(x)
            self.y = np.array(y)
            self.source_data["source"] = "5-digit"
            self.source_data["designation"] = designation
            return True

    def sourceAsXML(self, elemname):
        '''Saves airfoil source data in XML format
        
        Inputs:
            elemname: name of XML element to write
        Returns:
            elem: XML element
        '''

        elem = ET.Element(elemname)

        if self.source_data["source"] is None:
            return None

        ET.SubElement(elem, "source").text = self.source_data["source"]
        if self.source_data["source"] == "4-digit":
            ET.SubElement(elem, "camber").text = "{:.4f}".format(self.source_data["camber"])
            ET.SubElement(elem, "xCamber").text = "{:.4f}".format(self.source_data["xcamber"])
            ET.SubElement(elem, "thickness").text = "{:.4f}".format(self.source_data["thickness"])
        elif self.source_data["source"] == "5-digit":
            ET.SubElement(elem, "designation").text = self.source_data["designation"]
        elif self.source_data["source"] == "file":
            ET.SubElement(elem, "file").text = self.source_data["file"]

        return elem

    def fromXML(self, elem):
        '''Generates airfoil from source data in XML
        
        Inputs:
            elem: XML element containing airfoil source data
        Returns:
            retval: True on success, False on error
            errmsg: Description of error
        '''

        retval = True
        errmsg = "Success"
        source_data = {"source": None, "camber": None, "xcamber": None,
                       "thickness": None, "designation": None, "file": None}

        for child in elem:
            if child.tag == "source":
                source_data["source"] = child.text
            elif child.tag == "camber":
                source_data["camber"] = float(child.text)
            elif child.tag == "xCamber":
                source_data["xcamber"] = float(child.text)
            elif child.tag == "thickness":
                source_data["thickness"] = float(child.text)
            elif child.tag == "designation":
                source_data["designation"] = child.text
            elif child.tag == "file":
                source_data["file"] = child.text

        # Attempt to load airfoil
        if source_data["source"] == "4-digit":
            self.generate4Digit(source_data["camber"], source_data["xcamber"],
                                source_data["thickness"])
        elif source_data["source"] == "5-digit":
            retval = self.generate5Digit(source_data["designation"])
            if not retval:
                errmsg = "Invalid 5-digit seed airfoil designation."
        elif source_data["source"] == "file":
            check = self.readFromFile(source_data["file"])
            if check == 1:
                retval = False
                errmsg = "Failed to open {:s}.".format(source_data["file"])
            elif check == 2:
                retval = False
                errmsg = "Format error in airfoil file {:s}.".format(source_data["file"])

        return retval, errmsg

    def findLE(self):
        '''Finds leading edge as well as point location and where to split
        
        Returns:
            retval: True on success, False on error
            errmsg: message describing the error
        '''
        retval, errmsg = Airfoil.findLE(self)
        if not retval:
            return retval, errmsg

        # Determine leading edge index and where to add a point
        npt = self.numPoints()
        r1 = np.zeros((2))
        r2 = np.zeros((2))
        for i in range(npt-1):
            r1 = np.array([self.xle-self.x[i], self.yle-self.y[i]])
            dist1 = np.linalg.norm(r1)
            if dist1 != 0.:
                r1 /= dist1
            r2 = np.array([self.xle-self.x[i+1], self.yle-self.y[i+1]])
            dist2 = np.linalg.norm(r2)
            if dist2 != 0.:
                r2 /= dist2
            dot = np.dot(r1, r2)

            # Leading edge is between i and i+1 if any of these conditions occur
            if dist1 < 1e-12:
                self.leidx = i
                self.addpoint_loc = 0
                break
            elif dist2 < 1e-12:
                self.leidx = i
                self.addpoint_loc = 1
                break
            elif dot <= 0.:
                if dist1 < dist2:
                    self.leidx = i
                    self.addpoint_loc = 1
                else:
                    self.leidx = i+1
                    self.addpoint_loc = -1
                break

        # This shouldn't happen, but just in case...
        if self.leidx is None:
            retval = False
            errmsg = "Failed to determine leading edge point index. Please report this issue."

        return retval, errmsg

    def setOptimizationLimits(self, domaint, domainb):
        '''Sets optimization limits on top and bottom surfaces

        Inputs:
            domaint, domainb: upper and lower surface optimization domains

        Returns:
            retval: True on success, False on error
            errmsg: message describing the error
        '''
        errmsg = "Success"
        for domain in [domaint, domainb]:
            if domain[1] - domain[0] < 0.:
                errmsg = "Optimization domain must be positive."
                return False, errmsg
            elif domain[0] < 0.:
                errmsg = "Optimization lower limit cannot be less than 0."
                return False, errmsg
            elif domain[1] > 1.:
                errmsg = "Optimization uppwer limit cannot be greater than 1."
                return False, errmsg
        self.domaint = domaint
        self.domainb = domainb
        retval, errmsg = self.setDomainIdx()

        return retval, errmsg

    def split(self, symmetrical=False):
        '''Splits between upper and lower half
        
        Returns:
            retval: True on success, False on error
            errmsg: message describing the error'''

        # Determine number of points on upper and lower halves
        if self.addpoint_loc == 0:
            pointst = self.leidx+1
            pointsb = self.numPoints() - (self.leidx+1) + 1
            boundst = self.leidx-1
            boundsb = self.leidx+1
        elif self.addpoint_loc == -1:
            pointst = self.leidx+1
            pointsb = self.numPoints() - (self.leidx+1) + 2
            boundst = self.leidx-1
            boundsb = self.leidx
        else:
            pointst = self.leidx+2
            pointsb = self.numPoints() - (self.leidx+1) + 1
            boundst = self.leidx
            boundsb = self.leidx+1
        if symmetrical:
            pointsb = pointst

        # Split
        self.xt = np.zeros((pointst))
        self.yt = np.zeros((pointst))
        self.xt[0] = self.xle
        self.yt[0] = self.yle
        for i in range(pointst-1):
            self.xt[i+1] = self.x[boundst-i]
            self.yt[i+1] = self.y[boundst-i]

        self.xb = np.zeros((pointsb))
        self.yb = np.zeros((pointsb))
        self.xb[0] = self.xle
        self.yb[0] = self.yle
        if not symmetrical:
            for i in range(pointsb-1):
                self.xb[i+1] = self.x[boundsb+i]
                self.yb[i+1] = self.y[boundsb+i]
        else:
            for i in range(pointsb-1):
                self.xb[i+1] = self.xt[i+1]
                self.yb[i+1] = -self.yt[i+1]

        retval, errmsg = self.setDomainIdx()
        return retval, errmsg

    def setDomainIdx(self):
        '''Sets upper and lower optimization domain indices

        Returns:
            retval: True on success, False on error
            errmsg: message describing the error
        '''
        errmsg = "Success"

        # If airfoil hasn't been split yet, return as success. Will be set later when airfoil
        # is split.
        nptt = self.xt.shape[0] 
        nptb = self.xb.shape[0] 
        if (nptt == 0) and (nptb == 0):
            return True, errmsg
        elif (nptt == 0) or (nptb == 0):
            errmsg = "Airfoil split is not valid."
            return False, errmsg

        self.domainidxt = [0, 0]
        for i in range(nptt):
            xl = self.xt[self.domainidxt[0]]
            xr = self.xt[self.domainidxt[1]]
            if (self.xt[i] >= self.domaint[0]) and (xl < self.domaint[0]):
                self.domainidxt[0] = i
            if (self.xt[i] <= self.domaint[1]) and (xr < self.xt[i]):
                self.domainidxt[1] = i

        self.domainidxb = [0, 0]
        for i in range(nptb):
            xl = self.xb[self.domainidxb[0]]
            xr = self.xb[self.domainidxb[1]]
            if (self.xb[i] >= self.domainb[0]) and (xl < self.domainb[0]):
                self.domainidxb[0] = i
            if (self.xb[i] <= self.domainb[1]) and (xr < self.xb[i]):
                self.domainidxb[1] = i

        return True, errmsg

    def setupShapeFunctions(self, nshapest, nshapesb):
        '''Creates shape functions

        Inputs:
            nshapest, nshapesb: number of shape functions on top and bottom surfaces

        Returns:
            retval: True on success, False on error
            errmsg: Message describing the error
        '''
        errmsg = "Success"
        nptt = self.xt.shape[0]
        nptb = self.xb.shape[0]
        if (nptt == 0) or (nptb == 0):
            errmsg = "Airfoil must be split into upper and lower surfaces first."
            return False, errmsg

        if (self.domainidxt is None) or (self.domainidxb is None):
            errmsg = "Airfoil must have optimization limits set first."
            return False, errmsg

        self.shapest = []
        for i in range(nshapest):
            self.shapest.append(ShapeFunction(self.xt, self.domainidxt))
        self.shapesb = []
        for i in range(nshapesb):
            self.shapesb.append(ShapeFunction(self.xb, self.domainidxb))

        return True, errmsg

    def createPerturbedAirfoil(self, modest, modesb, t1fact, t2fact, symmetrical):
        '''Creates an airfoil using Hicks-Henne shape functions as perturbation

        Inputs:
            modest: design variables to define top surface shape functions
            modesb: design variables to define bottom surface shape functions
            t1fact: scaling factor for t1 design variable
            t2fact: scaling factor for t2 design variable
            symmetrical: whether the airfoil is symmetrical

        Returns:
            foil: perturbed airfoil
        '''

        # Account for scaling factors and create top surface
        yt = np.zeros((self.xt.shape[0]))
        for i, shape in enumerate(self.shapest):
            st = modest[3*i]
            t1 = modest[3*i+1]/t1fact
            t2 = modest[3*i+2]/t2fact
            self.shapest[i].createShape(st, t1, t2)
            yt += self.shapest[i].shape

        # Account for scaling factors and create bottom surface
        yb = np.zeros((self.xb.shape[0]))
        for i, shape in enumerate(self.shapesb):
            st = modesb[3*i]
            t1 = modesb[3*i+1]/t1fact
            t2 = modesb[3*i+2]/t2fact
            self.shapesb[i].createShape(st, t1, t2)
            yb += self.shapesb[i].shape

        # Create airfoil from top and bottom surfaces
        npt = self.xt.shape[0] + self.xb.shape[0] - 1
        foil = Airfoil(npt)
        foil.x = np.append(self.xt[::-1], self.xb[1:])
        foil.y = np.append(yt[::-1], yb[1:])

        return foil
