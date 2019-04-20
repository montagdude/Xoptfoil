import os
from PyQt5.QtWidgets import QDialog, QFileDialog
from settings import plotsettings
import plotsettings_ui

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
