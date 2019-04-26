from PyQt5.QtWidgets import QDialog

import flaphingepoint_ui

class FlapHingePointDialog(QDialog):
    def __init__(self):
        super(FlapHingePointDialog, self).__init__()
        self.ui = flaphingepoint_ui.Ui_Dialog()
        self.ui.setupUi(self)

    def populate(self, operatingpoints):
        self.ui.xflapBox.setValue(operatingpoints.xflap)
        self.ui.yflapBox.setValue(operatingpoints.yflap)
        self.ui.yflapSpecBox.setCurrentText(operatingpoints.yflapSpecification)

    def xflap(self):
        return self.ui.xflapBox.value()

    def yflap(self):
        return self.ui.yflapBox.value()

    def yflapSpecification(self):
        return self.ui.yflapSpecBox.currentText()
