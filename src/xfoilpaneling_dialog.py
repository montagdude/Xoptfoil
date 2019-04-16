from PyQt5.QtWidgets import QDialog
import xfoilpanelingsettings

class XfoilPanelingDialog(QDialog):
    def __init__(self):
        super(XfoilPanelingDialog, self).__init__()
        self.ui = xfoilpanelingsettings.Ui_Dialog()
        self.ui.setupUi(self)
