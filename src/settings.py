import os
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

def writePrettyXML(elem, f, indentlevel=0, indent='  '):
    '''Writes an XML element to file with nice indentation. xml.dom.minidom is used to do the
       writing, because etree.ElementTree dumps it all out on a single line without indentation.

    Inputs:
        elem: an ElementTree XML element
        f: an open file object for writing
        indentlevel: how many tab levels the entire element should be indented
        indent: indent character(s)

    Reference:
    https://stackoverflow.com/questions/17402323/use-xml-etree-elementtree-to-print-nicely-formatted-xml-files
    '''
    ugly = ET.tostring(elem)
    pretty = minidom.parseString(ugly).toprettyxml(indent=indent)
    prettylist = pretty.split('\n')
    for line in prettylist:
        # Skip XML header and blank lines
        if not line.startswith("<?xml version") and len(line) > 0:
            f.write(indentlevel*indent + line + "\n")

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
        self.autosaveFrequency = 100
        self.autosaveBasename = "optfoil"

    # Writes settings to open xml file
    def writeXML(self, f, indentlevel=0, indent='  '):
        root = ET.Element("OptimizationSettings")
        ET.SubElement(root, "searchType").text = self.searchType
        ET.SubElement(root, "globalSearch").text = self.globalSearch
        ET.SubElement(root, "localSearch").text = self.localSearch
        ET.SubElement(root, "shapeFunctions").text = self.shapeFunctions
        ET.SubElement(root, "nfunctionsTop").text = "{:d}".format(self.nfunctionsTop)
        ET.SubElement(root, "nfunctionsBot").text = "{:d}".format(self.nfunctionsBot)
        ET.SubElement(root, "initialPerturb").text = "{:.4f}".format(self.initialPerturb)
        ET.SubElement(root, "minBumpWidth").text = "{:.4f}".format(self.minBumpWidth)
        ET.SubElement(root, "autosaveFrequency").text = "{:d}".format(self.autosaveFrequency)
        ET.SubElement(root, "autosaveBasename").text = self.autosaveBasename
        writePrettyXML(root, f, indentlevel, indent)


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


class SimplexSettings():
    def __init__(self):
        self.maxIterations = 1000
        self.tolerance = 1e-06


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
simplexsettings = SimplexSettings()
xfoilsettings = XfoilSettings()
xfoilpanelingsettings = XfoilPanelingSettings()
plotsettings = PlotSettings()
