import os
import xml.etree.ElementTree as ET
import sys
from libxfoil import xfoil_geom_options_type

class Setting():
    def __init__(self, name=None, value=None, default=None, writeformat=None, datatype=None):
        self.name = name
        self.value = value
        self.default = default
        self.writeformat = writeformat
        self.datatype = datatype

        # Apply default value
        if (self.value is None) and (self.default is not None):
            self.value = self.default

class Settings():
    def __init__(self):
        self.settings = []

    def setting(self, name):
        for setting in self.settings:
            if setting.name == name:
                return setting
        return None

    def value(self, name):
        setting = self.setting(name)
        if setting is not None:
            return setting.value
        else:
            return None

    def addSetting(self, setting):
        if self.setting(setting.name) is None:
            self.settings.append(setting)

    def asXML(self, elemname):
        elem = ET.Element(elemname)
        for setting in self.settings:
            if setting.name is None:
                continue
            if setting.value is None:
                if setting.default is None:
                    text = "None"
                else:
                    text = setting.writeformat.format(setting.default)
            else:
                text = setting.writeformat.format(setting.value)
            ET.SubElement(elem, setting.name).text = text
        return elem

    def fromXML(self, elem):
        for child in elem:
            for setting in self.settings:
                if child.tag == setting.name:
                    if setting.datatype == bool:
                        setting.value = child.text == "True"
                    else:
                        try:
                            setting.value = setting.datatype(child.text)
                        except ValueError:
                            sys.stderr.write("Failed to read {:s}:{:s}: wrong data type.\n"\
                                             .format(elem.tag, setting.name))


class OptimizationSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="searchType", default="Global and local",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="globalSearch", default="Particle swarm",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="localSearch", default="Simplex", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="shapeFunctions", default="Hicks-Henne",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="nfunctionsTop", default=4, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="nfunctionsBot", default=4, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="initialPerturb", default=0.025, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="minBumpWidth", default=0.1, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="autosaveFrequency", default=100, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="autosaveBasename", default="optfoil", writeformat="{:s}",
                                datatype=str))


class InitializationSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="feasibleInit", default=True, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="feasibleLimit", default=50000., writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="maxAttempts", default=1000, writeformat="{:d}",
                                datatype=int))


class ParticleSwarmSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="population", default=40, writeformat="{:d}", datatype=int))
        self.addSetting(Setting(name="maxIterations", default=700, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="tolerance", default=1e-4, writeformat="{:.4e}",
                                datatype=float))
        self.addSetting(Setting(name="convergenceProfile", default="Exhaustive",
                                writeformat="{:s}", datatype=str))


class GeneticAlgorithmSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="population", default=80, writeformat="{:d}", datatype=int))
        self.addSetting(Setting(name="maxIterations", default=700, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="tolerance", default=1e-4, writeformat="{:.4e}",
                                datatype=float))
        self.addSetting(Setting(name="parentsSelection", default="Tournament",
                                writeformat="{:s}", datatype=str))
        self.addSetting(Setting(name="parentFraction", default=0.5, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="selectionPressure", default=8.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="tournamentFraction", default=0.025, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="crossoverFactor", default=0.5, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="mutantProbability", default=0.4, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="mutationRate", default=0.01, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="mutationFactor", default=0.2, writeformat="{:.4f}",
                                datatype=float))


class SimplexSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="maxIterations", default=1000, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="tolerance", default=1e-6, writeformat="{:.4e}",
                                datatype=float))


class XfoilSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="viscous", default=True, writeformat="{}", datatype=bool))
        self.addSetting(Setting(name="ncrit", default=9.0, writeformat="{:.4f}", datatype=float))
        self.addSetting(Setting(name="maxit", default=100, writeformat="{:d}", datatype=int))
        self.addSetting(Setting(name="xtript", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="xtripb", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="vaccel", default=0.01, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="fixUnconverged", default=True, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="reinitialize", default=True, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="silent", default=True, writeformat="{}", datatype=bool))
        self.silent = True


class XfoilPanelingSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="npan", default=160, writeformat="{:d}", datatype=int))
        self.addSetting(Setting(name="cvpar", default=1.0, writeformat="{:.4f}", datatype=float))
        self.addSetting(Setting(name="cterat", default=0.15, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="ctrrat", default=0.2, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="xsref1", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="xsref2", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="xpref1", default=1.0, writeformat="{:.4f}",
                                datatype=float))
        self.addSetting(Setting(name="xpref2", default=1.0, writeformat="{:.4f}",
                                datatype=float))

    def toXfoilGeomOpts(self):
        geom_opts = xfoil_geom_options_type()
        geom_opts.npan = self.value("npan")
        geom_opts.cvpar = self.value("cvpar")
        geom_opts.cterat = self.value("cterat")
        geom_opts.ctrrat = self.value("ctrrat")
        geom_opts.xsref1 = self.value("xsref1")
        geom_opts.xsref2 = self.value("xsref2")
        geom_opts.xpref1 = self.value("xpref1")
        geom_opts.xpref2 = self.value("xpref2")

        return geom_opts


class PlotSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
        self.addSetting(Setting(name="showSeedAirfoil", default=True, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="displayMode", default="Current design", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="designNumber", default=1, writeformat="{:d}",
                                datatype=int))
        self.addSetting(Setting(name="saveAnimationFrames", default=False, writeformat="{}",
                                datatype=bool))
        self.addSetting(Setting(name="frameDirectory", default=os.getcwd(), writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="framePrefix", default="optfoil", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="bgColor", default="#000000", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="fgColor", default="#ffffff", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="seedColor", default="#4b97ff", writeformat="{:s}",
                                datatype=str))
        self.addSetting(Setting(name="currentColor", default="#ff2424", writeformat="{:s}",
                                datatype=str))


# Instantiations
optimizationsettings = OptimizationSettings()
initializationsettings = InitializationSettings()
particleswarmsettings = ParticleSwarmSettings()
geneticalgorithmsettings = GeneticAlgorithmSettings()
simplexsettings = SimplexSettings()
xfoilsettings = XfoilSettings()
xfoilpanelingsettings = XfoilPanelingSettings()
plotsettings = PlotSettings()
