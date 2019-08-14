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
from constraints import constraints
from data import data
import methods
from settings import optimizationsettings

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
        self.ui.actionLoad_saved_settings.triggered.connect(self.loadSettings)
        self.ui.actionSave_settings.triggered.connect(self.saveSettings)
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

        self.ui.actionOptimize.triggered.connect(self.optimize)
        self.ui.actionPause.triggered.connect(self.pause)
        self.ui.actionStop.triggered.connect(self.stop)

    # Loads saved settings from XML file
    def loadSettings(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load saved settings",
                                                         self._currpath, "XML files (*.xml)")
        if fname == "":
            return
        else:
            self._currpath = os.path.dirname(fname)

        ret, errmsg = methods.read_settings(fname)
        if ret != 0:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
        else:
            # Airfoil and/or plot settings may have changed
            if data.seed_airfoil.numPoints() > 0:
                self.ui.mplwidget.plotAirfoils()
            else:
                self.ui.mplwidget.setupAxes()
                self.ui.mplwidget.draw()

            # Enable optimization if possible
            if methods.can_optimize():
                self.setOptimizationEnabled(True)

    # Saves settings to XML file
    def saveSettings(self):
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save settings", self._currpath,
                                                         "XML files (*.xml)")
        if fname == "":
            return
        else:
            self._currpath = os.path.dirname(fname)

        if not methods.save_settings(fname):
            QtWidgets.QMessageBox.critical(self, "Error",
                                           "Could not open {:s} for writing.".format(fname))

    # Loads seed airfoil from file
    def loadSeed(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load seed airfoil",
                                                         self._currpath)
        if fname == "":
            return
        else:
            self._currpath = os.path.dirname(fname)

        retval, errmsg = data.readSeedAirfoil(fname)
        if not retval:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
        else:
            self.ui.mplwidget.plotAirfoils()

            # Enable optimization if possible
            if methods.can_optimize():
                self.setOptimizationEnabled(True)

    # Generates a NACA seed airfoil
    def generateNACA(self):
        if self.naca_dialog.exec():
            retval = True
            if self.naca_dialog.airfoilType() == "4 digit":
                data.generate4DigitAirfoil(self.naca_dialog.camber(), self.naca_dialog.xcamber(),
                                           self.naca_dialog.thickness())
            elif self.naca_dialog.airfoilType() == "5 digit":
                retval = data.generate5DigitAirfoil(self.naca_dialog.designation())
            if retval:
                self.ui.mplwidget.plotAirfoils()

                # Enable optimization if possible
                if methods.can_optimize():
                    self.setOptimizationEnabled(True)
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

            # Enable optimization if possible
            if methods.can_optimize():
                self.setOptimizationEnabled(True)

    def setConstraints(self):
        dialog = ConstraintsDialog()
        if dialog.exec():
            dialog.saveConstraints()

    def setOptimizationEnabled(self, enabled):
        self.ui.optimizeButton.setEnabled(enabled)
        self.ui.actionOptimize.setEnabled(enabled)

    def optimize(self):
        # Enable/disable actions
        self.ui.pauseButton.setEnabled(True)
        self.ui.actionPause.setEnabled(True)
        self.ui.stopButton.setEnabled(True)
        self.ui.actionStop.setEnabled(True)
        self.setOptimizationEnabled(False)

        # Find LE and transform seed airfoil
        retval, errmsg = data.processSeedAirfoil()
        if not retval:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
            return

        # Update plot to reflect seed airfoil transformation
        self.ui.mplwidget.plotAirfoils()

        # Set optimization limits
        domaint = [optimizationsettings.value("topLeftLimit"),
                   optimizationsettings.value("topRightLimit")]
        domainb = [optimizationsettings.value("bottomLeftLimit"),
                   optimizationsettings.value("bottomRightLimit")]
        retval, errmsg = data.seed_airfoil.setOptimizationLimits(domaint, domainb)
        if not retval:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
            return

        # Split seed airfoil
        retval, errmsg = data.seed_airfoil.split(constraints.value("symmetrical"))
        if not retval:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
            return

        # Setup shape functions
        nshapest = optimizationsettings.value("nfunctionsTop")
        nshapesb = optimizationsettings.value("nfunctionsBot")
        shapetype = optimizationsettings.value("shapeFunctions")
        retval, errmsg = data.seed_airfoil.setupShapeFunctions(nshapest, nshapesb, shapetype)
        if not retval:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
            return

        # Check that seed airfoil passes constraints
        _, failures = data.seed_airfoil.checkGeometry(constraints, data.xoffset, data.foilscale)
        if len(failures) > 0:
            msg = "Seed airfoil violates the following constraints:\n"
            for failure in failures:
                msg += "    " + failure + "\n"
            msg += "\nContinue optimizing anyway?"
            reply = QtWidgets.QMessageBox.question(self, "Seed airfoil constraint violations",
                                                   msg, QtWidgets.QMessageBox.Yes |
                                                        QtWidgets.QMessageBox.No,
                                                        QtWidgets.QMessageBox.No)
            if reply == QtWidgets.QMessageBox.No:
                return

    def pause(self):
        # Enable/disable actions
        self.ui.pauseButton.setEnabled(False)
        self.ui.actionPause.setEnabled(False)
        self.ui.stopButton.setEnabled(True)
        self.ui.actionStop.setEnabled(True)
        self.setOptimizationEnabled(True)

    def stop(self):
        # Enable/disable actions
        self.ui.pauseButton.setEnabled(False)
        self.ui.actionPause.setEnabled(False)
        self.ui.stopButton.setEnabled(False)
        self.ui.actionStop.setEnabled(False)
        self.setOptimizationEnabled(True)


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    application = XoptfoilMainWindow()
    application.show()
    sys.exit(app.exec())
