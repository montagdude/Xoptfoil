from PyQt5.QtWidgets import QDialog
import xfoilpanelingsettings_ui
from settings import xfoilpanelingsettings

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
