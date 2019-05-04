from PyQt5.QtWidgets import QDialog, QFileDialog, QColorDialog, QMessageBox
from PyQt5.QtGui import QRegExpValidator, QDoubleValidator, QColor
from PyQt5.QtCore import QRegExp

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

        # RegExp validator for basename line edit (any number of alphanumeric characters)
        regexp = QRegExp("^[a-zA-Z0-9]*")
        regexp_validator = QRegExpValidator(regexp)
        self.ui.basenameEdit.setValidator(regexp_validator)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.searchTypeBox.setCurrentText(optimizationsettings.value("searchType"))
        self.ui.globalSearchBox.setCurrentText(optimizationsettings.value("globalSearch"))
        self.ui.localSearchBox.setCurrentText(optimizationsettings.value("localSearch"))
        self.ui.shapeFunctionsBox.setCurrentText(optimizationsettings.value("shapeFunctions"))
        self.ui.nfunctionsTopBox.setValue(optimizationsettings.value("nfunctionsTop"))
        self.ui.nfunctionsBotBox.setValue(optimizationsettings.value("nfunctionsBot"))
        self.ui.initialPerturbBox.setValue(optimizationsettings.value("initialPerturb"))
        self.ui.minBumpWidthBox.setValue(optimizationsettings.value("minBumpWidth"))
        self.ui.autosaveBox.setValue(optimizationsettings.value("autosaveFrequency"))
        self.ui.basenameEdit.setText(optimizationsettings.value("autosaveBasename"))

    def saveSettings(self):
        optimizationsettings.setting("searchType").value = self.ui.searchTypeBox.currentText()
        optimizationsettings.setting("globalSearch").value = self.ui.globalSearchBox.currentText()
        optimizationsettings.setting("localSearch").value = self.ui.localSearchBox.currentText()
        optimizationsettings.setting("shapeFunctions").value = \
                             self.ui.shapeFunctionsBox.currentText()
        optimizationsettings.setting("nfunctionsTop").value = self.ui.nfunctionsTopBox.value()
        optimizationsettings.setting("nfunctionsBot").value = self.ui.nfunctionsBotBox.value()
        optimizationsettings.setting("initialPerturb").value = self.ui.initialPerturbBox.value()
        optimizationsettings.setting("minBumpWidth").value = self.ui.minBumpWidthBox.value()
        optimizationsettings.setting("autosaveFrequency").value = self.ui.autosaveBox.value()
        optimizationsettings.setting("autosaveBasename").value = self.ui.basenameEdit.text()


class InitializationSettingsDialog(QDialog):
    def __init__(self):
        super(InitializationSettingsDialog, self).__init__()
        self.ui = initializationsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.feasibleInitBox.setChecked(initializationsettings.value("feasibleInit"))
        self.ui.feasibleLimitBox.setValue(initializationsettings.value("feasibleLimit"))
        self.ui.maxAttemptsBox.setValue(initializationsettings.value("maxAttempts"))

    def saveSettings(self):
        initializationsettings.setting("feasibleInit").value = self.ui.feasibleInitBox.isChecked()
        initializationsettings.setting("feasibleLimit").value = self.ui.feasibleLimitBox.value()
        initializationsettings.setting("maxAttempts").value = self.ui.maxAttemptsBox.value()


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
        self.ui.populationBox.setValue(particleswarmsettings.value("population"))
        self.ui.maxIterationsBox.setValue(particleswarmsettings.value("maxIterations"))
        self.ui.toleranceEdit.setText("{:.4e}".format(particleswarmsettings.value("tolerance")))
        self.ui.convergenceProfileBox.setCurrentText(
                                      particleswarmsettings.value("convergenceProfile"))

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
        particleswarmsettings.setting("population").value = self.ui.populationBox.value()
        particleswarmsettings.setting("maxIterations").value = self.ui.maxIterationsBox.value()
        particleswarmsettings.setting("tolerance").value = float(self.ui.toleranceEdit.text())
        particleswarmsettings.setting("convergenceProfile").value = \
                              self.ui.convergenceProfileBox.currentText()


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
        self.ui.populationBox.setValue(geneticalgorithmsettings.value("population"))
        self.ui.maxIterationsBox.setValue(geneticalgorithmsettings.value("maxIterations"))
        self.ui.toleranceEdit.setText("{:.4e}".format(
                              geneticalgorithmsettings.value("tolerance")))
        self.ui.parentsSelectionBox.setCurrentText(
                                    geneticalgorithmsettings.value("parentsSelection"))
        self.ui.parentFractionBox.setValue(geneticalgorithmsettings.value("parentFraction"))
        self.ui.selectionPressureBox.setValue(geneticalgorithmsettings.value("selectionPressure"))
        self.ui.tournamentFractionBox.setValue(
                                      geneticalgorithmsettings.value("tournamentFraction"))
        self.ui.crossoverFactorBox.setValue(geneticalgorithmsettings.value("crossoverFactor"))
        self.ui.mutantProbabilityBox.setValue(geneticalgorithmsettings.value("mutantProbability"))
        self.ui.mutationRateBox.setValue(geneticalgorithmsettings.value("mutationRate"))
        self.ui.mutationFactorBox.setValue(geneticalgorithmsettings.value("mutationFactor"))

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
        geneticalgorithmsettings.setting("population").value = self.ui.populationBox.value()
        geneticalgorithmsettings.setting("maxIterations").value = self.ui.maxIterationsBox.value()
        geneticalgorithmsettings.setting("tolerance").value = float(self.ui.toleranceEdit.text())
        geneticalgorithmsettings.setting("parentsSelection").value = \
                                 self.ui.parentsSelectionBox.currentText()
        geneticalgorithmsettings.setting("parentFraction").value = \
                                 self.ui.parentFractionBox.value()
        geneticalgorithmsettings.setting("selectionPressure").value = \
                                 self.ui.selectionPressureBox.value()
        geneticalgorithmsettings.setting("tournamentFraction").value = \
                                 self.ui.tournamentFractionBox.value()
        geneticalgorithmsettings.setting("crossoverFactor").value = \
                                 self.ui.crossoverFactorBox.value()
        geneticalgorithmsettings.setting("mutantProbability").value = \
                                 self.ui.mutantProbabilityBox.value()
        geneticalgorithmsettings.setting("mutationRate").value = self.ui.mutationRateBox.value()
        geneticalgorithmsettings.setting("mutationFactor").value = \
                                 self.ui.mutationFactorBox.value()


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
        self.ui.maxIterationsBox.setValue(simplexsettings.value("maxIterations"))
        self.ui.toleranceEdit.setText("{:.4e}".format(simplexsettings.value("tolerance")))

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
        simplexsettings.setting("maxIterations").value = self.ui.maxIterationsBox.value()
        simplexsettings.setting("tolerance").value = float(self.ui.toleranceEdit.text())


class XfoilSettingsDialog(QDialog):
    def __init__(self):
        super(XfoilSettingsDialog, self).__init__()
        self.ui = xfoilsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.viscousBox.setChecked(xfoilsettings.value("viscous"))
        self.ui.ncritBox.setValue(xfoilsettings.value("ncrit"))
        self.ui.maxitBox.setValue(xfoilsettings.value("maxit"))
        self.ui.xtriptBox.setValue(xfoilsettings.value("xtript"))
        self.ui.xtripbBox.setValue(xfoilsettings.value("xtripb"))
        self.ui.vaccelBox.setValue(xfoilsettings.value("vaccel"))
        self.ui.fixUnconvergedBox.setChecked(xfoilsettings.value("fixUnconverged"))
        self.ui.reinitializeBox.setChecked(xfoilsettings.value("reinitialize"))
        self.ui.silentBox.setChecked(xfoilsettings.value("silent"))

    def saveSettings(self):
        xfoilsettings.setting("viscous").value = self.ui.viscousBox.isChecked()
        xfoilsettings.setting("ncrit").value = self.ui.ncritBox.value()
        xfoilsettings.setting("maxit").value = self.ui.maxitBox.value()
        xfoilsettings.setting("xtript").value = self.ui.xtriptBox.value()
        xfoilsettings.setting("xtripb").value = self.ui.xtripbBox.value()
        xfoilsettings.setting("vaccel").value = self.ui.vaccelBox.value()
        xfoilsettings.setting("fixUnconverged").value = self.ui.fixUnconvergedBox.isChecked()
        xfoilsettings.setting("reinitialize").value = self.ui.reinitializeBox.isChecked()
        xfoilsettings.setting("silent").value = self.ui.silentBox.isChecked()


class XfoilPanelingSettingsDialog(QDialog):
    def __init__(self):
        super(XfoilPanelingSettingsDialog, self).__init__()
        self.ui = xfoilpanelingsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        self.ui.npanBox.setValue(xfoilpanelingsettings.value("npan"))
        self.ui.cvparBox.setValue(xfoilpanelingsettings.value("cvpar"))
        self.ui.cteratBox.setValue(xfoilpanelingsettings.value("cterat"))
        self.ui.ctrratBox.setValue(xfoilpanelingsettings.value("ctrrat"))
        self.ui.xsref1Box.setValue(xfoilpanelingsettings.value("xsref1"))
        self.ui.xsref2Box.setValue(xfoilpanelingsettings.value("xsref2"))
        self.ui.xpref1Box.setValue(xfoilpanelingsettings.value("xpref1"))
        self.ui.xpref2Box.setValue(xfoilpanelingsettings.value("xpref2"))

    def saveSettings(self):
        xfoilpanelingsettings.setting("npan").value = self.ui.npanBox.value()
        xfoilpanelingsettings.setting("cvpar").value = self.ui.cvparBox.value()
        xfoilpanelingsettings.setting("cterat").value = self.ui.cteratBox.value()
        xfoilpanelingsettings.setting("ctrrat").value = self.ui.ctrratBox.value()
        xfoilpanelingsettings.setting("xsref1").value = self.ui.xsref1Box.value()
        xfoilpanelingsettings.setting("xsref2").value = self.ui.xsref2Box.value()
        xfoilpanelingsettings.setting("xpref1").value = self.ui.xpref1Box.value()
        xfoilpanelingsettings.setting("xpref2").value = self.ui.xpref2Box.value()


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
        self.bgColor = plotsettings.value("bgColor")
        self.fgColor = plotsettings.value("fgColor")
        self.seedColor = plotsettings.value("seedColor")
        self.currentColor = plotsettings.value("currentColor")

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
                         "Animation frame directory", plotsettings.value("frameDirectory"))
        self.ui.frameDirectoryEdit.setText(frameDirectory)

    def setBgColor(self):
        newcolor = QColorDialog.getColor(QColor(plotsettings.value("bgColor")))
        if newcolor.isValid():
            self.bgColor = newcolor.name()
            self.ui.bgColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                                self.bgColor))

    def setFgColor(self):
        newcolor = QColorDialog.getColor(QColor(plotsettings.value("fgColor")))
        if newcolor.isValid():
            self.fgColor = newcolor.name()
            self.ui.fgColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                                self.fgColor))

    def setSeedColor(self):
        newcolor = QColorDialog.getColor(QColor(plotsettings.value("seedColor")))
        if newcolor.isValid():
            self.seedColor = newcolor.name()
            self.ui.seedColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                                  self.seedColor))

    def setCurrentColor(self):
        newcolor = QColorDialog.getColor(QColor(plotsettings.value("currentColor")))
        if newcolor.isValid():
            self.currentColor = newcolor.name()
            self.ui.currentColorButton.setStyleSheet("background-color:{:s}; border:0px;".format(
                                                     self.currentColor))

    def populate(self):
        self.ui.showSeedBox.setChecked(plotsettings.value("showSeedAirfoil"))
        self.ui.displayModeBox.setCurrentText(plotsettings.value("displayMode"))
        self.ui.designNumberBox.setValue(plotsettings.value("designNumber"))
        self.ui.saveFrameBox.setChecked(plotsettings.value("saveAnimationFrames"))
        self.ui.frameDirectoryEdit.setText(plotsettings.value("frameDirectory"))
        self.ui.framePrefixEdit.setText(plotsettings.value("framePrefix"))

    def saveSettings(self):
        plotsettings.setting("showSeedAirfoil").value = self.ui.showSeedBox.isChecked()
        plotsettings.setting("displayMode").value = self.ui.displayModeBox.currentText()
        plotsettings.setting("designNumber").value = self.ui.designNumberBox.value()
        plotsettings.setting("saveAnimationFrames").value = self.ui.saveFrameBox.isChecked()
        plotsettings.setting("frameDirectory").value = self.ui.frameDirectoryEdit.text()
        plotsettings.setting("framePrefix").value = self.ui.framePrefixEdit.text()
        plotsettings.setting("bgColor").value = self.bgColor
        plotsettings.setting("fgColor").value = self.fgColor
        plotsettings.setting("seedColor").value = self.seedColor
        plotsettings.setting("currentColor").value = self.currentColor
