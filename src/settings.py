import os

class OptimizationSettings():
    def __init__(self):
        self.searchType = "Global and local"
        self.globalSearch = "Particle swarm"
        self.localSearch = "Simplex"
        self.shapeFunctions = "Hicks-Henne"
        self.nfunctionsTop = 4
        self.nfunctionsBot = 4
        self.initialPerturb = 0.025
        self.minBumpWidth = 0.1


class InitializationSettings():
    def __init__(self):
        self.feasibleInit = True
        self.feasibleLimit = 50000.
        self.maxAttempts = 1000


class ParticleSwarmSettings():
    def __init__(self):
        self.population = 40
        self.maxIterations = 700
        self.tolerance = 1e-4
        self.convergenceProfile = "Exhaustive"


class XfoilSettings():
    def __init__(self):
        self.viscous = True
        self.ncrit = 9.0
        self.maxit = 100
        self.xtript = 1.0
        self.xtripb = 1.0
        self.vaccel = 0.01
        self.fixUnconverged = True
        self.reinitialize = True
        self.silent = True


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


class PlotSettings():
    def __init__(self):
        self.showSeedAirfoil = True
        self.displayMode = "Current design"
        self.designNumber = 1
        self.saveAnimationFrames = False
        self.frameDirectory = os.getcwd()
        self.framePrefix = "optfoil"


# Instantiations
optimizationsettings = OptimizationSettings()
initializationsettings = InitializationSettings()
particleswarmsettings = ParticleSwarmSettings()
xfoilsettings = XfoilSettings()
xfoilpanelingsettings = XfoilPanelingSettings()
plotsettings = PlotSettings()
