from PyQt5.QtWidgets import QDialog, QFileDialog
import optimizationsettings_ui
import initializationsettings_ui
import particleswarmsettings_ui
import xfoilsettings_ui
import plotsettings_ui
import xfoilpanelingsettings_ui
from settings import (optimizationsettings, initializationsettings, particleswarmsettings,
                      xfoilsettings, xfoilpanelingsettings, plotsettings)

class OptimizationSettingsDialog(QDialog):
    def __init__(self):
        super(OptimizationSettingsDialog, self).__init__()
        self.ui = optimizationsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.searchTypeBox.setCurrentText(optimizationsettings.searchType)
        self.ui.globalSearchBox.setCurrentText(optimizationsettings.globalSearch)
        self.ui.localSearchBox.setCurrentText(optimizationsettings.localSearch)
        self.ui.shapeFunctionsBox.setCurrentText(optimizationsettings.shapeFunctions)
        self.ui.nfunctionsTopBox.setValue(optimizationsettings.nfunctionsTop)
        self.ui.nfunctionsBotBox.setValue(optimizationsettings.nfunctionsBot)
        self.ui.initialPerturbBox.setValue(optimizationsettings.initialPerturb)
        self.ui.minBumpWidthBox.setValue(optimizationsettings.minBumpWidth)

    def saveSettings(self):
        optimizationsettings.searchType = self.ui.searchTypeBox.currentText()
        optimizationsettings.globalSearch = self.ui.globalSearchBox.currentText()
        optimizationsettings.localSearch = self.ui.localSearchBox.currentText()
        optimizationsettings.shapeFunctions = self.ui.shapeFunctionsBox.currentText()
        optimizationsettings.nfunctionsTop = self.ui.nfunctionsTopBox.value()
        optimizationsettings.nfunctionsBot = self.ui.nfunctionsBotBox.value()
        optimizationsettings.initialPerturb = self.ui.initialPerturbBox.value()
        optimizationsettings.minBumpWidth = self.ui.minBumpWidthBox.value()


class InitializationSettingsDialog(QDialog):
    def __init__(self):
        super(InitializationSettingsDialog, self).__init__()
        self.ui = initializationsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.feasibleInitBox.setChecked(initializationsettings.feasibleInit)
        self.ui.feasibleLimitBox.setValue(initializationsettings.feasibleLimit)
        self.ui.maxAttemptsBox.setValue(initializationsettings.maxAttempts)

    def saveSettings(self):
        initializationsettings.feasibleInit = self.ui.feasibleInitBox.isChecked()
        initializationsettings.feasibleLimit = self.ui.feasibleLimitBox.value()
        initializationsettings.maxAttempts = self.ui.maxAttemptsBox.value()


class ParticleSwarmSettingsDialog(QDialog):
    def __init__(self):
        super(ParticleSwarmSettingsDialog, self).__init__()
        self.ui = particleswarmsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.populationBox.setValue(particleswarmsettings.population)
        self.ui.maxIterationsBox.setValue(particleswarmsettings.maxIterations)
        self.ui.toleranceBox.setValue(particleswarmsettings.tolerance)
        self.ui.convergenceProfileBox.setCurrentText(particleswarmsettings.convergenceProfile)

    def saveSettings(self):
        particleswarmsettings.population = self.ui.populationBox.value()
        particleswarmsettings.maxIterations = self.ui.maxIterationsBox.value()
        particleswarmsettings.tolerance = self.ui.toleranceBox.value()
        particleswarmsettings.convergenceProfile = self.ui.convergenceProfileBox.currentText()


class XfoilSettingsDialog(QDialog):
    def __init__(self):
        super(XfoilSettingsDialog, self).__init__()
        self.ui = xfoilsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.viscousBox.setChecked(xfoilsettings.viscous)
        self.ui.ncritBox.setValue(xfoilsettings.ncrit)
        self.ui.maxitBox.setValue(xfoilsettings.maxit)
        self.ui.xtriptBox.setValue(xfoilsettings.xtript)
        self.ui.xtripbBox.setValue(xfoilsettings.xtripb)
        self.ui.vaccelBox.setValue(xfoilsettings.vaccel)
        self.ui.fixUnconvergedBox.setChecked(xfoilsettings.fixUnconverged)
        self.ui.reinitializeBox.setChecked(xfoilsettings.reinitialize)
        self.ui.silentBox.setChecked(xfoilsettings.silent)

    def saveSettings(self):
        xfoilsettings.viscous = self.ui.viscousBox.isChecked()
        xfoilsettings.ncrit = self.ui.ncritBox.value()
        xfoilsettings.maxit = self.ui.maxitBox.value()
        xfoilsettings.xtript = self.ui.xtriptBox.value()
        xfoilsettings.xtripb = self.ui.xtripbBox.value()
        xfoilsettings.vaccel = self.ui.vaccelBox.value()
        xfoilsettings.fixUnconverged = self.ui.fixUnconvergedBox.isChecked()
        xfoilsettings.reinitialize = self.ui.reinitializeBox.isChecked()
        xfoilsettings.silent = self.ui.silentBox.isChecked()


class XfoilPanelingSettingsDialog(QDialog):
    def __init__(self):
        super(XfoilPanelingSettingsDialog, self).__init__()
        self.ui = xfoilpanelingsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.npanBox.setValue(xfoilpanelingsettings.npan)
        self.ui.cvparBox.setValue(xfoilpanelingsettings.cvpar)
        self.ui.cteratBox.setValue(xfoilpanelingsettings.cterat)
        self.ui.ctrratBox.setValue(xfoilpanelingsettings.ctrrat)
        self.ui.xsref1Box.setValue(xfoilpanelingsettings.xsref1)
        self.ui.xsref2Box.setValue(xfoilpanelingsettings.xsref2)
        self.ui.xpref1Box.setValue(xfoilpanelingsettings.xpref1)
        self.ui.xpref2Box.setValue(xfoilpanelingsettings.xpref2)

    def saveSettings(self):
        xfoilpanelingsettings.npan = self.ui.npanBox.value()
        xfoilpanelingsettings.cvpar = self.ui.cvparBox.value()
        xfoilpanelingsettings.cterat = self.ui.cteratBox.value()
        xfoilpanelingsettings.ctrrat = self.ui.ctrratBox.value()
        xfoilpanelingsettings.xsref1 = self.ui.xsref1Box.value()
        xfoilpanelingsettings.xsref2 = self.ui.xsref2Box.value()
        xfoilpanelingsettings.xpref1 = self.ui.xpref1Box.value()
        xfoilpanelingsettings.xpref2 = self.ui.xpref2Box.value()


class PlotSettingsDialog(QDialog):
    def __init__(self):
        super(PlotSettingsDialog, self).__init__()
        self.ui = plotsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Signals/slots
        self.ui.displayModeBox.currentIndexChanged[str].connect(self.setItemsEnabled)
        self.ui.browseButton.clicked.connect(self.setAnimationDirectory)

        # Populate items
        self.populate()

    # Enable/disable widgets based on display mode
    def setItemsEnabled(self, choice):
        single_design_widgets = [self.ui.designNumberLabel, self.ui.designNumberBox]
        animate_widgets = [self.ui.saveFrameBox, self.ui.frameDirectoryLabel,
                           self.ui.frameDirectoryEdit, self.ui.browseButton,
                           self.ui.framePrefixLabel, self.ui.framePrefixEdit]
        if choice == "Current design":
            for widget in single_design_widgets:
                widget.setEnabled(False)
            for widget in animate_widgets:
                widget.setEnabled(False)
        elif choice == "Single design":
            for widget in single_design_widgets:
                widget.setEnabled(True)
            for widget in animate_widgets:
                widget.setEnabled(False)
        elif choice == "Animate all designs":
            for widget in single_design_widgets:
                widget.setEnabled(False)
            for widget in animate_widgets:
                widget.setEnabled(True)

    def setAnimationDirectory(self):
        frameDirectory = QFileDialog.getExistingDirectory(self,
                         "Animation frame directory", plotsettings.frameDirectory)
        self.ui.frameDirectoryEdit.setText(frameDirectory)

    def populate(self):
        self.ui.showSeedBox.setChecked(plotsettings.showSeedAirfoil)
        self.ui.displayModeBox.setCurrentText(plotsettings.displayMode)
        self.ui.designNumberBox.setValue(plotsettings.designNumber)
        self.ui.saveFrameBox.setChecked(plotsettings.saveAnimationFrames)
        self.ui.frameDirectoryEdit.setText(plotsettings.frameDirectory)
        self.ui.framePrefixEdit.setText(plotsettings.framePrefix)

    def saveSettings(self):
        plotsettings.showSeedAirfoil = self.ui.showSeedBox.isChecked()
        plotsettings.displayMode = self.ui.displayModeBox.currentText()
        plotsettings.designNumber = self.ui.designNumberBox.value()
        plotsettings.saveAnimationFrames = self.ui.saveFrameBox.isChecked()
        plotsettings.frameDirectory = self.ui.frameDirectoryEdit.text()
        plotsettings.framePrefix = self.ui.framePrefixEdit.text()
