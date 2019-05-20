import sys
import numpy as np
import xml.etree.ElementTree as ET
import libxfoil_wrap as lxf
from libxfoil import xfoil_data_group

class Airfoil:

    npointside = 100

    def __init__(self):
        self.x = np.zeros((0))
        self.y = np.zeros((0))
        self.xdg = xfoil_data_group()
        lxf.xfoil_init(self.xdg)
        self.source_data = {"source": None, "camber": None, "xcamber": None,
                            "thickness": None, "designation": None, "file": None}
        self.xle = None
        self.yle = None

    def numPoints(self):
        return self.x.shape[0]

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
        if len(x) > 1:
            self.x = np.array(x)
            self.y = np.array(y)
        else:
            f.close()
            errmsg = "Less than 2 vertices in coordinates file."
            sys.stderr.write(errmsg+"\n")
            return False, errmsg
        f.close()

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
        errmsg = "Success"
        if self.numPoints() < 3:
            errmsg = "Airfoil must have at least 3 points."
            return False, errmsg
        s, xs, ys = lxf.xfoil_spline_coordinates(self.x, self.y, self.numPoints())
        _, self.xle, self.yle = lxf.xfoil_lefind(self.x, self.y, s, xs, ys, self.numPoints())

        return True, errmsg

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

        return np.amax(self.x) - self.xle, True, errmsg

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
