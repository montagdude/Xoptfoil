from PyQt5.QtWidgets import QDialog
import xfoilsettings_ui

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
