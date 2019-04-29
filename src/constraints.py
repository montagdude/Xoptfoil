import sys

class Constraints():
    def __init__(self):
        self.minThickness = 6.0
        self.maxThickness = 1000.0
        self.minCamber = -10.0
        self.maxCamber = 10.0
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
        self.addThickX = []
        self.addThickMin = []
        self.addThickMax = []

    def numMomentConstraints(self):
        return len(self.momentConstraintType)

    def addMomentConstraint(constraint_type, min_moment=None):
        self.momentConstraintType.append(constraint_type)
        if constraint_type == "Specify":
            if min_moment is None:
                sys.stderr.write("min_moment must be provided for specified constraint type.\n")
                sys.exit(1)
            self.minMoment.append(min_moment)
        else:
            self.minMoment.append(0.0)

    def numThicknessConstraints(self):
        return len(self.addThickX)

    def addThicknessConstraint(x, minthick, maxthick):
        self.addThickX.append(x)
        self.addThickMin.append(minthick)
        self.addThickMax.append(maxthick)

constraints = Constraints()
