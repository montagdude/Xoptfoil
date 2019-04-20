import os

class PlotSettings():
    def __init__(self):
        self.showSeedAirfoil = True
        self.displayMode = "Current design"
        self.designNumber = 1
        self.saveAnimationFrames = False
        self.frameDirectory = os.getcwd()
        self.framePrefix = "optfoil"

class XfoilPanelingSettings():
    def __init__(self):
        self.npan = 160
        self.cvpar = 1.0
        self.cterat = 0.15
        self.ctrrat = 0.2
        self.xsref1 = 1.0
        self.xsref2 = 1.0
        self.xpref1 = 1.0
        self.xpref2 = 1.0

# Instantiations
plotsettings = PlotSettings()
xfoilpanelingsettings = XfoilPanelingSettings()
