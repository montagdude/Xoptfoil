from PyQt5.QtWidgets import QDialog, QMessageBox
from PyQt5.QtGui import QDoubleValidator

import operatingpoint_ui
from operatingpoints import OperatingPoint
from settings import xfoilsettings

class OperatingPointDialog(QDialog):
    def __init__(self):
        super(OperatingPointDialog, self).__init__()
        self.ui = operatingpoint_ui.Ui_Dialog()
        self.ui.setupUi(self)

        self.updateCondition(self.ui.specConditionBox.currentText())

        # Double validator for Re line edit
        dbl_validator = QDoubleValidator(bottom=0.)

        # Signals and slots
        self.ui.specConditionBox.currentIndexChanged[str].connect(self.updateCondition)
        self.ui.flapBehaviorBox.currentIndexChanged[str].connect(self.setFlapVisibility)
        self.ui.ncritBehaviorBox.currentIndexChanged[str].connect(self.setNcritVisibility)
        self.ui.tripBehaviorBox.currentIndexChanged[str].connect(self.setTripVisibility)

    def updateCondition(self, choice):
        self.ui.conditionLabel.setText(choice)
        if choice == "Cl":
            self.ui.conditionBox.setSingleStep(0.1)
        elif choice == "AoA":
            self.ui.conditionBox.setSingleStep(1.0)

    def setFlapVisibility(self, choice):
        flap_widgets = [self.ui.flapDeflectionLabel, self.ui.flapDeflectionBox]
        if choice == "Specified deflection":
            for widget in flap_widgets:
                widget.setEnabled(True)
        else:
            for widget in flap_widgets:
                widget.setEnabled(False)

    def setNcritVisibility(self, choice):
        ncrit_widgets = [self.ui.ncritLabel, self.ui.ncritBox]
        if choice == "Use Xfoil settings":
            for widget in ncrit_widgets:
                widget.setEnabled(False)
        else:
            for widget in ncrit_widgets:
                widget.setEnabled(True)

    def setTripVisibility(self, choice):
        trip_widgets = [self.ui.xtriptLabel, self.ui.xtriptBox, self.ui.xtripbLabel,
                        self.ui.xtripbBox]
        if choice == "Use Xfoil settings":
            for widget in trip_widgets:
                widget.setEnabled(False)
        else:
            for widget in trip_widgets:
                widget.setEnabled(True)

    # Sets data in the form from a given operating point
    def fromOperatingPoint(self, point):
        self.ui.optimizationGoalBox.setCurrentText(point.value("optimizationGoal"))
        self.ui.specConditionBox.setCurrentText(point.value("specCondition"))
        self.ui.conditionBox.setValue(point.value("condition"))
        self.ui.reynoldsEdit.setText("{:.4e}".format(point.value("reynolds")))
        self.ui.machBox.setValue(point.value("mach"))
        self.ui.flapBehaviorBox.setCurrentText(point.value("flapBehavior"))
        self.ui.flapDeflectionBox.setValue(point.value("flapDeflection"))
        self.ui.ncritBehaviorBox.setCurrentText(point.value("ncritBehavior"))
        if point.value("ncritBehavior") == "Use Xfoil settings":
            self.ui.ncritBox.setValue(xfoilsettings.value("ncrit"))
        else:
            self.ui.ncritBox.setValue(point.value("ncrit"))
        self.ui.tripBehaviorBox.setCurrentText(point.value("tripBehavior"))
        if point.value("tripBehavior") == "Use Xfoil settings":
            self.ui.xtriptBox.setValue(xfoilsettings.value("xtript"))
            self.ui.xtripbBox.setValue(xfoilsettings.value("xtripb"))
        else:
            self.ui.xtriptBox.setValue(point.value("xtript"))
            self.ui.xtripbBox.setValue(point.value("xtripb"))
        self.ui.weightingBox.setValue(point.value("weighting"))

    # Returns data in the form of an operating point
    def operatingPoint(self):
        point = OperatingPoint()
        point.setting("optimizationGoal").value = self.ui.optimizationGoalBox.currentText()
        point.setting("specCondition").value = self.ui.specConditionBox.currentText()
        point.setting("condition").value = self.ui.conditionBox.value()
        point.setting("reynolds").value = float(self.ui.reynoldsEdit.text())
        point.setting("mach").value = self.ui.machBox.value()
        point.setting("flapBehavior").value = self.ui.flapBehaviorBox.currentText()
        point.setting("flapDeflection").value = self.ui.flapDeflectionBox.value()
        point.setting("ncritBehavior").value = self.ui.ncritBehaviorBox.currentText()
        if point.value("ncritBehavior") == "Use Xfoil settings":
            point.setting("ncrit").value = xfoilsettings.value("ncrit")
        else:
            point.setting("ncrit").value = self.ui.ncritBox.value()
        point.setting("tripBehavior").value = self.ui.tripBehaviorBox.currentText()
        if point.value("tripBehavior") == "Use Xfoil settings":
            point.setting("xtript").value = xfoilsettings.value("xtript")
            point.setting("xtripb").value = xfoilsettings.value("xtripb")
        else:
            point.setting("xtript").value = self.ui.xtriptBox.value()
            point.setting("xtripb").value = self.ui.xtripbBox.value()
        point.setting("weighting").value = self.ui.weightingBox.value()

        return point

    def accept(self):
        if self.validateSettings():
            self.setResult(QDialog.Accepted)
            self.hide()
        else:
            self.setResult(QDialog.Rejected)

    def validateSettings(self):
        try:
            check = float(self.ui.reynoldsEdit.text())
        except ValueError:
            QMessageBox.critical(self, "Error", "Cannot convet {:s} to float"\
                                 .format(self.ui.reynoldsEdit.text()))
            return False
        
        # Check for bad combinations of operating conditions and optimization goals
        specCondition = self.ui.specConditionBox.currentText()
        condition = self.ui.conditionBox.value()
        optimizationGoal = self.ui.optimizationGoalBox.currentText()
        if (specCondition == "Cl") and (optimizationGoal == "Maximize lift"):
            QMessageBox.critical(self, "Error", "Cannot maximize lift when Cl is specified.")
            return False
        elif (condition <= 0.0) and (specCondition == "Cl"):
            if (optimizationGoal == "Maximize lift") or \
               (optimizationGoal == "Maximize glide slope"):
                   QMessageBox.critical(self, "Error", "Cl must be greather 0 to use the "
                                        + "{:s} optimization goal.".format(optimizationGoal))
                   return False
        return True
