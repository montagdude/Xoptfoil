from PyQt5.QtWidgets import QDialog
import optimizationsettings_ui

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
