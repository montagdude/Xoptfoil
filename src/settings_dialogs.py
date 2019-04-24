from PyQt5.QtWidgets import QDialog, QFileDialog, QColorDialog, QMessageBox
from PyQt5.QtGui import QDoubleValidator, QColor

import optimizationsettings_ui
import initializationsettings_ui
import particleswarmsettings_ui
import geneticalgorithmsettings_ui
import simplexsettings_ui
import xfoilsettings_ui
import plotsettings_ui
import xfoilpanelingsettings_ui
from settings import (optimizationsettings, initializationsettings, particleswarmsettings,
                      geneticalgorithmsettings, simplexsettings, xfoilsettings,
                      xfoilpanelingsettings, plotsettings)

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

        # Double validator for tolerance line edit
        dbl_validator = QDoubleValidator(bottom=0.)
        self.ui.toleranceEdit.setValidator(dbl_validator)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.populationBox.setValue(particleswarmsettings.population)
        self.ui.maxIterationsBox.setValue(particleswarmsettings.maxIterations)
        self.ui.toleranceEdit.setText("{:.4e}".format(particleswarmsettings.tolerance))
        self.ui.convergenceProfileBox.setCurrentText(particleswarmsettings.convergenceProfile)

    def accept(self):
        if self.validateSettings():
            self.setResult(QDialog.Accepted)
            self.hide()
        else:
            self.setResult(QDialog.Rejected)

    def validateSettings(self):
        try:
            check = float(self.ui.toleranceEdit.text())
        except ValueError:
            QMessageBox.critical(self, "Error", "Cannot convert {:s} to float."\
                                 .format(self.ui.toleranceEdit.text()))
            return False
        return True

    def saveSettings(self):
        particleswarmsettings.population = self.ui.populationBox.value()
        particleswarmsettings.maxIterations = self.ui.maxIterationsBox.value()
        particleswarmsettings.tolerance = float(self.ui.toleranceEdit.text())
        particleswarmsettings.convergenceProfile = self.ui.convergenceProfileBox.currentText()


class GeneticAlgorithmSettingsDialog(QDialog):
    def __init__(self):
        super(GeneticAlgorithmSettingsDialog, self).__init__()
        self.ui = geneticalgorithmsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Double validator for tolerance line edit
        dbl_validator = QDoubleValidator(bottom=0.)
        self.ui.toleranceEdit.setValidator(dbl_validator)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.populationBox.setValue(geneticalgorithmsettings.population)
        self.ui.maxIterationsBox.setValue(geneticalgorithmsettings.maxIterations)
        self.ui.toleranceEdit.setText("{:.4e}".format(geneticalgorithmsettings.tolerance))
        self.ui.parentsSelectionBox.setCurrentText(geneticalgorithmsettings.parentsSelection)
        self.ui.parentFractionBox.setValue(geneticalgorithmsettings.parentFraction)
        self.ui.selectionPressureBox.setValue(geneticalgorithmsettings.selectionPressure)
        self.ui.tournamentFractionBox.setValue(geneticalgorithmsettings.tournamentFraction)
        self.ui.crossoverFactorBox.setValue(geneticalgorithmsettings.crossoverFactor)
        self.ui.mutantProbabilityBox.setValue(geneticalgorithmsettings.mutantProbability)
        self.ui.mutationRateBox.setValue(geneticalgorithmsettings.mutationRate)
        self.ui.mutationFactorBox.setValue(geneticalgorithmsettings.mutationFactor)

    def accept(self):
        if self.validateSettings():
            self.setResult(QDialog.Accepted)
            self.hide()
        else:
            self.setResult(QDialog.Rejected)

    def validateSettings(self):
        try:
            check = float(self.ui.toleranceEdit.text())
        except ValueError:
            QMessageBox.critical(self, "Error", "Cannot convert {:s} to float."\
                                 .format(self.ui.toleranceEdit.text()))
            return False
        return True

    def saveSettings(self):
        geneticalgorithmsettings.population = self.ui.populationBox.value()
        geneticalgorithmsettings.maxIterations = self.ui.maxIterationsBox.value()
        geneticalgorithmsettings.tolerance = float(self.ui.toleranceEdit.text())
        geneticalgorithmsettings.parentsSelection = self.ui.parentsSelectionBox.currentText()
        geneticalgorithmsettings.parentFraction = self.ui.parentFractionBox.value()
        geneticalgorithmsettings.selectionPressure = self.ui.selectionPressureBox.value()
        geneticalgorithmsettings.tournamentFraction = self.ui.tournamentFractionBox.value()
        geneticalgorithmsettings.crossoverFactor = self.ui.crossoverFactorBox.value()
        geneticalgorithmsettings.mutantProbability = self.ui.mutantProbabilityBox.value()
        geneticalgorithmsettings.mutationRate = self.ui.mutationRateBox.value()
        geneticalgorithmsettings.mutationFactor = self.ui.mutationFactorBox.value()


class SimplexSettingsDialog(QDialog):
    def __init__(self):
        super(SimplexSettingsDialog, self).__init__()
        self.ui = simplexsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Double validator for tolerance line edit
        dbl_validator = QDoubleValidator(bottom=0.)
        self.ui.toleranceEdit.setValidator(dbl_validator)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.maxIterationsBox.setValue(simplexsettings.maxIterations)
        self.ui.toleranceEdit.setText("{:.4e}".format(simplexsettings.tolerance))

    def accept(self):
        if self.validateSettings():
            self.setResult(QDialog.Accepted)
            self.hide()
        else:
            self.setResult(QDialog.Rejected)

    def validateSettings(self):
        try:
            check = float(self.ui.toleranceEdit.text())
        except ValueError:
            QMessageBox.critical(self, "Error", "Cannot convert {:s} to float."\
                                 .format(self.ui.toleranceEdit.text()))
            return False
        return True

    def saveSettings(self):
        simplexsettings.maxIterations = self.ui.maxIterationsBox.value()
        simplexsettings.tolerance = float(self.ui.toleranceEdit.text())


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
        self.ui.bgColorButton.clicked.connect(self.setBgColor)
        self.ui.fgColorButton.clicked.connect(self.setFgColor)
        self.ui.seedColorButton.clicked.connect(self.setSeedColor)
        self.ui.currentColorButton.clicked.connect(self.setCurrentColor)

        # Colors
        self.bgColor = plotsettings.bgColor
        self.fgColor = plotsettings.fgColor
        self.seedColor = plotsettings.seedColor
        self.currentColor = plotsettings.currentColor

        # Button colors
        self.ui.bgColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                            self.bgColor))
        self.ui.fgColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                            self.fgColor))
        self.ui.seedColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                              self.seedColor))
        self.ui.currentColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                                  self.currentColor))

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

    def setBgColor(self):
        self.bgColor = QColorDialog.getColor(QColor(plotsettings.bgColor)).name()
        self.ui.bgColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                            self.bgColor))

    def setFgColor(self):
        self.fgColor = QColorDialog.getColor(QColor(plotsettings.fgColor)).name()
        self.ui.fgColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                            self.fgColor))

    def setSeedColor(self):
        self.seedColor = QColorDialog.getColor(QColor(plotsettings.seedColor)).name()
        self.ui.seedColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                              self.seedColor))

    def setCurrentColor(self):
        self.currentColor = QColorDialog.getColor(QColor(plotsettings.currentColor)).name()
        self.ui.currentColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                                 self.currentColor))

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
        plotsettings.bgColor = self.bgColor
        plotsettings.fgColor = self.fgColor
        plotsettings.seedColor = self.seedColor
        plotsettings.currentColor = self.currentColor
