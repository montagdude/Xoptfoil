from PyQt5.QtWidgets import QDialog
import xfoilpanelingsettings_ui

class XfoilPanelingDialog(QDialog):
    def __init__(self):
        super(XfoilPanelingDialog, self).__init__()
        self.ui = xfoilpanelingsettings_ui.Ui_Dialog()
        self.ui.setupUi(self)
