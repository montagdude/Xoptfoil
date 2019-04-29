#!/usr/bin/python3

import sys
import os
from PyQt5 import QtWidgets

# FIXME: this will have to be set at install time. Also doesn't work if not
# run in the src directory, so obviously this is only temporary.
installdir = os.path.join(os.getcwd(), '..')
sys.path.append(os.path.join(installdir, 'ui'))

import mainwindow_ui
from nacagenerator_dialog import NacaGeneratorDialog
from settings_dialogs import (OptimizationSettingsDialog, InitializationSettingsDialog,
                              ParticleSwarmSettingsDialog, GeneticAlgorithmSettingsDialog,
                              SimplexSettingsDialog, XfoilSettingsDialog,
                              XfoilPanelingSettingsDialog, PlotSettingsDialog)
from operatingpoints_dialog import OperatingPointsDialog
from constraints_dialog import ConstraintsDialog
from data import data

class XoptfoilMainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(XoptfoilMainWindow, self).__init__()
        self._currpath = os.getcwd()

        self.ui = mainwindow_ui.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle("Xoptfoil")

        # NACA generator dialog. Make this persistent so settings are kept across uses.
        self.naca_dialog = NacaGeneratorDialog()

        # Connect toolbar to canvas
        self.ui.mpltoolbar.setCanvas(self.ui.mplwidget)

        # Signals and slots
        self.ui.actionLoad_seed_airfoil.triggered.connect(self.loadSeed)
        self.ui.actionGenerate_NACA_airfoil.triggered.connect(self.generateNACA)
        self.ui.actionQuit.triggered.connect(self.close)

        self.ui.actionOptimization_settings.triggered.connect(self.showOptimizationSettings)
        self.ui.actionInitialization_settings.triggered.connect(self.showInitializationSettings)
        self.ui.actionParticleswarm_settings.triggered.connect(self.showParticleSwarmSettings)
        self.ui.actionGenetic_algorithm_settings.triggered.connect(
                self.showGeneticAlgorithmSettings)
        self.ui.actionSimplex_settings.triggered.connect(self.showSimplexSettings)
        self.ui.actionXfoil_settings.triggered.connect(self.showXfoilSettings)
        self.ui.actionXfoil_paneling_settings.triggered.connect(self.showXfoilPanelingSettings)
        self.ui.actionPlot_settings.triggered.connect(self.showPlotSettings)

        self.ui.actionSet_Operating_Points.triggered.connect(self.setOperatingPoints)

        self.ui.actionSet_constraints.triggered.connect(self.setConstraints)

    # Loads seed airfoil from file
    def loadSeed(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load seed airfoil",
                                                         self._currpath)
        if fname == "":
            return
        else:
            self._currpath = os.path.dirname(fname)

        retval = data.readSeedAirfoil(fname)
        errmsg = None
        if retval == 1:
            errmsg  = "Unable to read airfoil file {:s}: I/O error.".format(fname)
        elif retval == 2:
            errmsg  = "Unable to read airfoil file {:s}:".format(fname)
            errmsg += " Not a valid airfoil file."
        if errmsg is not None:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
        self.ui.mplwidget.plotAirfoils()

    # Generates a NACA seed airfoil
    def generateNACA(self):
        if self.naca_dialog.exec():
            if self.naca_dialog.airfoilType() == "4 digit":
                retval = data.generate4DigitAirfoil(self.naca_dialog.camber(),
                                                    self.naca_dialog.xcamber(),
                                                    self.naca_dialog.thickness())
            elif self.naca_dialog.airfoilType() == "5 digit":
                retval = data.generate5DigitAirfoil(self.naca_dialog.designation())
            if retval:
                self.ui.mplwidget.plotAirfoils()
            else:
                QtWidgets.QMessageBox.critical(self, "Error", "Unsupported 5-digit designation.")

    def showOptimizationSettings(self):
        dialog = OptimizationSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showInitializationSettings(self):
        dialog = InitializationSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showParticleSwarmSettings(self):
        dialog = ParticleSwarmSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showGeneticAlgorithmSettings(self):
        dialog = GeneticAlgorithmSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showSimplexSettings(self):
        dialog = SimplexSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showXfoilSettings(self):
        dialog = XfoilSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showXfoilPanelingSettings(self):
        dialog = XfoilPanelingSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()

    def showPlotSettings(self):
        dialog = PlotSettingsDialog()
        if dialog.exec():
            dialog.saveSettings()
            # FIXME: do this judiciously and plot the thing that is actually selected
            if len(data.seed_airfoil.x) > 0:
                self.ui.mplwidget.plotAirfoils()
            else:
                self.ui.mplwidget.setupAxes()
                self.ui.mplwidget.draw()

    def setOperatingPoints(self):
        dialog = OperatingPointsDialog()
        if dialog.exec():
            dialog.saveOperatingPoints()

    def setConstraints(self):
        dialog = ConstraintsDialog()
        if dialog.exec():
            dialog.saveConstraints()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    application = XoptfoilMainWindow()
    application.show()
    sys.exit(app.exec())
