import sys
import xml.etree.ElementTree as ET

from settings import Settings, Setting

class Constraints(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="minThickness", default=0.06, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="maxThickness", default=100.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="minCamber", default=-0.1, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="maxCamber", default=0.1, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="minTEAngle", default=3.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="checkCurvature", default=True, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="maxReverseTop", default=0, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="maxReverseBot", default=1, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="curveThreshold", default=0.075, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="symmetrical", default=False, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="minFlapAngle", default=-20.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="maxFlapAngle", default=20.0, writeformat="{:.4f}",
                                datatype=float))
        self.momentConstraintType = []
        self.minMoment = []
        self.cpConstraintType = []
        self.minCp = []
        self.addThickX = []
        self.addThickMin = []
        self.addThickMax = []

    def numMomentConstraints(self):
        return len(self.momentConstraintType)

    def addMomentConstraint(self, constraint_type, min_moment=None):
        self.momentConstraintType.append(constraint_type)
        if constraint_type == "Specify":
            if min_moment is None:
                sys.stderr.write("min_moment must be provided for specified constraint type.\n")
                sys.exit(1)
            self.minMoment.append(min_moment)
        else:
            self.minMoment.append(0.0)

    def clearMomentConstraints(self):
        self.momentConstraintType = []
        self.minMoment = []

    def numCpConstraints(self):
        return len(self.cpConstraintType)

    def addCpConstraint(self, constraint_type, min_cp=None):
        self.cpConstraintType.append(constraint_type)
        if constraint_type == "Specify":
            if min_cp is None:
                sys.stderr.write("min_cp must be provided for specified constraint type.\n")
                sys.exit(1)
            self.minCp.append(min_cp)
        else:
            self.minCp.append(0.0)

    def clearCpConstraints(self):
        self.cpConstrainteType = []
        self.minCp = []

    def numThicknessConstraints(self):
        return len(self.addThickX)

    def addThicknessConstraint(self, x, minthick, maxthick):
        self.addThickX.append(x)
        self.addThickMin.append(minthick)
        self.addThickMax.append(maxthick)

    def clearThicknessConstraints(self):
        self.addThickX = []
        self.addThickMin = []
        self.addThickMax = []

    def asXML(self, elemname):
        elem = Settings.asXML(self, elemname)

        # Moment constraints
        nmoment = self.numMomentConstraints()
        for i in range(nmoment):
            subelem = ET.Element("MomentConstraint")
            ET.SubElement(subelem, "momentConstraintType").text = self.momentConstraintType[i]
            if self.momentConstraintType[i] == "Specify":
                ET.SubElement(subelem, "minMoment").text = "{:.4f}".format(self.minMoment[i])
            elem.append(subelem)

        # Cp constraints
        ncp = self.numCpConstraints()
        for i in range(ncp):
            subelem = ET.Element("CpConstraint")
            ET.SubElement(subelem, "cpConstraintType").text = self.cpConstraintType[i]
            if self.cpConstraintType[i] == "Specify":
                ET.SubElement(subelem, "minCp").text = "{:.4f}".format(self.minCp[i])
            elem.append(subelem)

        # Additional thickness constraints
        nthick = self.numThicknessConstraints()
        for i in range(nthick):
            subelem = ET.Element("ThicknessConstraint")
            ET.SubElement(subelem, "xThickness").text = "{:.4f}".format(self.addThickX[i])
            ET.SubElement(subelem, "minThickness").text = "{:.4f}".format(self.addThickMin[i])
            ET.SubElement(subelem, "maxThickness").text = "{:.4f}".format(self.addThickMax[i])
            elem.append(subelem)

        return elem

    def fromXML(self, elem):
        Settings.fromXML(self, elem)
        self.clearMomentConstraints()
        self.clearCpConstraints()
        self.clearThicknessConstraints()

        for child in elem:
            # Read moment constraints
            if child.tag == "MomentConstraint":
                for grandchild in child:
                    min_moment = None
                    if grandchild.tag == "momentConstraintType":
                        constraint_type = grandchild.text
                    elif grandchild.tag == "minMoment":
                        min_moment = float(grandchild.text)
                self.addMomentConstraint(constraint_type, min_moment)

            # Read Cp constraints
            elif child.tag == "CpConstraint":
                for grandchild in child:
                    min_cp = None
                    if grandchild.tag == "cpConstraintType":
                        constraint_type = grandchild.text
                    elif grandchild.tag == "minCp":
                        min_cp = float(grandchild.text)
                self.addCpConstraint(constraint_type, min_cp)

            # Read thickness constraints
            elif child.tag == "ThicknessConstraint":
                for grandchild in child:
                    if grandchild.tag == "xThickness":
                        x = float(grandchild.text)
                    elif grandchild.tag == "minThickness":
                        minthick = float(grandchild.text)
                    elif grandchild.tag == "maxThickness":
                        maxthick = float(grandchild.text)
                # Do some basic validation and add it
                if x < 0. or x > 1.:
                    sys.stderr.write("Warning: xThickness must be >= 0 and <= 1. Skipping.\n")
                elif minthick >= maxthick:
                    sys.stderr.write("Warning: minThickness must be < maxThickness. Skipping.\n")
                else:
                    self.addThicknessConstraint(x, minthick, maxthick)


constraints = Constraints()
