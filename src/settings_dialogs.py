from PyQt5.QtWidgets import QDialog, QFileDialog
import optimizationsettings_ui
import xfoilsettings_ui
import plotsettings_ui
import xfoilpanelingsettings_ui
from settings import xfoilpanelingsettings, plotsettings

class OptimizationSettingsDialog(QDialog):
    def __init__(self):
        super(OptimizationSettingsDialog, self).__init__()
        self.ui = optimizationsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        pass

    def saveSettings(self):
        pass


class XfoilSettingsDialog(QDialog):
    def __init__(self):
        super(XfoilSettingsDialog, self).__init__()
        self.ui = xfoilsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Populate items
        self.populate()

    def populate(self):
        pass

    def saveSettings(self):
        pass


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
