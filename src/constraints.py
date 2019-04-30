import sys

class Constraints():
    def __init__(self):
        self.minThickness = 0.06
        self.maxThickness = 100.0
        self.minCamber = -0.1
        self.maxCamber = 0.1
        self.minTEAngle = 3.0
        self.checkCurvature = True
        self.maxReverseTop = 0
        self.maxReverseBot = 1
        self.curveThreshold = 0.075
        self.symmetrical = False
        self.minFlapAngle = -20.
        self.maxFlapAngle = 20.
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

constraints = Constraints()
