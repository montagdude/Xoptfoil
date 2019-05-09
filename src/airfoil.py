import sys
import numpy as np
import xml.etree.ElementTree as ET
import libxfoil_wrap as lxf

class Airfoil:
    def __init__(self):
        self.x = []
        self.y = []
        self.source_data = {"source": None, "camber": None, "xcamber": None,
                            "thickness": None, "designation": None, "file": None}


    def readFromFile(self, fname):
        """Reads airfoil from file.

        Inputs
            fname: file name or path
        Returns
            0 on success
            1 on file IO error
            2 on format error
        """
        
        print("Reading airfoil file {:s}".format(fname))
        try:
            f = open(fname)
        except IOError:
            sys.stderr.write("Unable to read {:s}\n.".format(fname))
            return 1

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
                    sys.stderr.write("Not enough data in line.\n")
                    return 2
                try:
                    x.append(float(splitline[0]))
                    y.append(float(splitline[1]))
                except ValueError:
                    f.close()
                    sys.stderr.write("Unable to convert to float.\n")
                    return 2
            linenum += 1
        if len(x) > 1:
            self.x = np.array(x)
            self.y = np.array(y)
        else:
            f.close()
            sys.stderr.write("Less than 2 vertices in coordinates file.\n")
            return 2
        f.close()

        self.source_data["source"] = "file"
        self.source_data["file"] = fname
        return 0

    def generate4Digit(self, camber, xcamber, thickness, npointside):
        x, y, _ = lxf.naca_4_digit(camber, xcamber, thickness, npointside)
        self.x = np.array(x)
        self.y = np.array(y)
        self.source_data["source"] = "4-digit"
        self.source_data["camber"] = camber
        self.source_data["xcamber"] = xcamber
        self.source_data["thickness"] = thickness

    def generate5Digit(self, designation, npointside):
        x, y, _, stat = lxf.naca_5_digit(designation, npointside)
        if stat != 0:
            return False
        else:
            self.x = np.array(x)
            self.y = np.array(y)
            self.source_data["source"] = "5-digit"
            self.source_data["designation"] = designation
            return True

    def sourceAsXML(self, elemname):
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
