import xml.etree.ElementTree as ET

from settings import Settings, Setting

class OperatingPoint(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="optimizationGoal", default="Minimize drag",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="specCondition", default="Cl", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="condition", default=0.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="reynolds", default=1e+05, writeformat="{:.4e}",
                                datatype=float))
        self.addSetting(Setting(name="mach", default=0.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="flapBehavior", default="No deflection", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="flapDeflection", default=0.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="ncritBehavior", default="Use Xfoil settings",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="ncrit", default=9.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="tripBehavior", default="Use Xfoil settings",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="xtript", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="xtripb", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="weighting", default=1.0, writeformat="{:.4f}",
                                datatype=float))


class OperatingPoints(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="xflap", default=0.7, writeformat="{:.4f}", datatype=float))
        self.addSetting(Setting(name="yflapSpecification", default="y / c",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="yflap", default=0.0, writeformat="{:.4f}", datatype=float))
        self.points = []

    def numPoints(self):
        return len(self.points)

    def point(self, idx):
        if idx < self.numPoints():
            return self.points[idx]
        else:
            return None

    def setPoint(self, idx, point):
        if idx < self.numPoints():
            self.points[idx] = point

    def addPoint(self, point):
        self.points.append(point)

    def deletePoint(self, idx):
        if idx < self.numPoints():
            self.points.pop(idx)

    def deleteAllPoints(self):
        self.points = []

    def asXML(self, elemname):
        elem = Settings.asXML(self, elemname)
        ET.SubElement(elem, "numPoints").text = "{:d}".format(self.numPoints())
        for point in self.points:
            elem.append(point.asXML("OperatingPoint"))

        return elem


operatingpoints = OperatingPoints()
