from PyQt5.QtWidgets import QDialog, QMessageBox

import nacagenerator_ui

class NacaGeneratorDialog(QDialog):
    def __init__(self):
        super(NacaGeneratorDialog, self).__init__()
        self.ui = nacagenerator_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Signals/slots
        self.ui.typeBox.currentIndexChanged[str].connect(self.setItemsEnabled)

    # Enable/disable widgets based on airfoil type
    def setItemsEnabled(self, choice):
        four_digit_widgets = [self.ui.camberLabel, self.ui.camberBox, self.ui.xcamberLabel,
                              self.ui.xcamberBox, self.ui.thicknessLabel, self.ui.thicknessBox]
        five_digit_widgets = [self.ui.designationLabel, self.ui.designationEdit]
        if choice == "4 digit":
            for widget in four_digit_widgets:
                widget.setEnabled(True)
            for widget in five_digit_widgets:
                widget.setEnabled(False)
        elif choice == "5 digit":
            for widget in four_digit_widgets:
                widget.setEnabled(False)
            for widget in five_digit_widgets:
                widget.setEnabled(True)

    def airfoilType(self):
        return self.ui.typeBox.currentText()

    def camber(self):
        return self.ui.camberBox.value()

    def xcamber(self):
        return self.ui.xcamberBox.value()

    def thickness(self):
        return self.ui.thicknessBox.value()

    def designation(self):
        return self.ui.designationEdit.text()

    def accept(self):
        if self.validateSettings():
            self.setResult(QDialog.Accepted)
            self.hide()
        else:
            self.setResult(QDialog.Rejected)

    def validateSettings(self):
        if self.ui.typeBox.currentText() == "4 digit":
            return True
        des = self.ui.designationEdit.text()
        valid_firstthree = ['210', '220', '230', '240', '250']
        if not des[0:3] in valid_firstthree:
            QMessageBox.critical(self, "Error", "First three digits must be " +
                                 "210, 220, 230, 240, or 250.")
            return False
        thickmsg = "Last two digits must be a thickness from 00 to 99."
        try:
            thick = float(des[3:])
        except ValueError:
            QMessageBox.critical(self, "Error", thickmsg)
            return False
        if thick < 0.0 or thick > 99.0:
            QMessageBox.critical(self, "Error", thickmsg)
            return False
        return True
