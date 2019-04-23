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


class GeneticAlgorithmSettings():
    def __init__(self):
        self.population = 80
        self.maxIterations = 700
        self.tolerance = 1e-04
        self.parentsSelection = "Tournament"
        self.parentFraction = 0.5
        self.selectionPressure = 8.0
        self.tournamentFraction = 0.025
        self.crossoverFactor = 0.5
        self.mutantProbability = 0.4
        self.mutationRate = 0.01
        self.mutationFactor = 0.2


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
        self.bgColor = '#000000'
        self.fgColor = '#ffffff'
        self.seedColor = '#4b97ff'
        self.currentColor = '#ff2424'


# Instantiations
optimizationsettings = OptimizationSettings()
initializationsettings = InitializationSettings()
particleswarmsettings = ParticleSwarmSettings()
geneticalgorithmsettings = GeneticAlgorithmSettings()
xfoilsettings = XfoilSettings()
xfoilpanelingsettings = XfoilPanelingSettings()
plotsettings = PlotSettings()
