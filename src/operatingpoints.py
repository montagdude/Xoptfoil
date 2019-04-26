class OperatingPoint():
    def __init__(self):
        self.specCondition = "Cl"
        self.condition = 0.0
        self.reynolds = 1e+05
        self.mach = 0.0
        self.flapBehavior = "No deflection"
        self.flapDeflection = 0.0
        self.ncritBehavior = "Use Xfoil settings"
        self.ncrit = 9.0
        self.tripBehavior = "Use Xfoil settings"
        self.xtript = 1.0
        self.xtripb = 1.0
        self.weighting = 1.0

class OperatingPoints():
    def __init__(self):
        self.xflap = 0.7
        self.yflapSpecification = "y / c"
        self.yflap = 0.0
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

operatingpoints = OperatingPoints()
