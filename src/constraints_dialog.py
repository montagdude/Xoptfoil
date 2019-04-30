from PyQt5.QtWidgets import (QDialog, QMessageBox, QHeaderView, QTableWidgetItem, QComboBox,
                             QDoubleSpinBox)

import constraints_ui
from constraints import constraints
from operatingpoints import operatingpoints

class ConstraintsDialog(QDialog):
    def __init__(self):
        super(ConstraintsDialog, self).__init__()
        self.ui = constraints_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Resize behavior for headers
        self.ui.momentTable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.ui.cpTable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.ui.thicknessTable.horizontalHeader().setSectionResizeMode(
                                                  QHeaderView.ResizeToContents)

        # Signals and slots
        self.ui.checkCurvatureBox.stateChanged[int].connect(self.setItemsEnabled)
        self.ui.addButton.clicked.connect(self.addThicknessConstraint)
        self.ui.deleteButton.clicked.connect(self.deleteThicknessConstraint)

        # Sets up moment constraints table
        self.momenttype_boxes = []
        self.minmoment_boxes = []
        self.setupMomentTable()

        # Sets up cp constraints table
        self.cptype_boxes = []
        self.mincp_boxes = []
        self.setupCpTable()

        # Additional thickness constraints widgets
        self.addthickx_boxes = []
        self.addthickmin_boxes = []
        self.addthickmax_boxes = []

        # Populate items
        self.populate()

    # Populates items
    def populate(self):
        self.ui.minThicknessBox.setValue(constraints.minThickness*100.)
        self.ui.maxThicknessBox.setValue(constraints.maxThickness*100.)
        self.ui.minCamberBox.setValue(constraints.minCamber*100.)
        self.ui.maxCamberBox.setValue(constraints.maxCamber*100.)
        self.ui.minTEAngleBox.setValue(constraints.minTEAngle)
        self.ui.checkCurvatureBox.setChecked(constraints.checkCurvature)
        self.ui.maxReverseTopBox.setValue(constraints.maxReverseTop)
        self.ui.maxReverseBottomBox.setValue(constraints.maxReverseBot)
        self.ui.curveThresholdBox.setValue(constraints.curveThreshold)
        self.ui.symmetricalBox.setChecked(constraints.symmetrical)
        self.ui.minFlapAngleBox.setValue(constraints.minFlapAngle)
        self.ui.maxFlapAngleBox.setValue(constraints.maxFlapAngle)

        # Min moment boxes. Make consistent with number of operating points, then set values.
        numpoints = operatingpoints.numPoints()
        if len(self.momenttype_boxes) > numpoints:
            self.momenttype_boxes = self.momenttype_boxes[0:numpoints]
            self.minmoment_boxes = self.minmoment_boxes[0:numpoints]
        elif len(self.momenttype_boxes) < numpoints:
            nadd = numpoints - len(self.momenttype_boxes)
            for i in range(nadd):
                self.appendMomentBoxes()
        for i in range(min(constraints.numMomentConstraints(),numpoints)):
            self.momenttype_boxes[i].setCurrentText(constraints.momentConstraintType[i])
            self.minmoment_boxes[i].setValue(constraints.minMoment[i])

        # Min Cp boxes. Treat same way as moment boxes.
        if len(self.cptype_boxes) > numpoints:
            self.cptype_boxes = self.cptype_boxes[0:numpoints]
            self.mincp_boxes = self.mincp_boxes[0:numpoints]
        elif len(self.cptype_boxes) < numpoints:
            nadd = numpoints - len(self.cptype_boxes)
            for i in range(nadd):
                self.appendCpBoxes()
        for i in range(min(constraints.numCpConstraints(),numpoints)):
            self.cptype_boxes[i].setCurrentText(constraints.cpConstraintType[i])
            self.mincp_boxes[i].setValue(constraints.minCp[i])

        # Thickness constraints
        for i in range(constraints.numThicknessConstraints()):
            self.addThicknessConstraint()
            self.addthickx_boxes[i].setValue(constraints.addThickX[i])
            self.addthickmin_boxes[i].setValue(constraints.addThickMin[i]*100.)
            self.addthickmax_boxes[i].setValue(constraints.addThickMax[i]*100.)

    # Enables/disables items
    def setItemsEnabled(self, choice):
        if choice == 2:
            self.ui.maxReverseTopBox.setEnabled(True)
            self.ui.maxReverseBottomBox.setEnabled(True)
            self.ui.curveThresholdBox.setEnabled(True)
        elif choice == 0:
            self.ui.maxReverseTopBox.setEnabled(False)
            self.ui.maxReverseBottomBox.setEnabled(False)
            self.ui.curveThresholdBox.setEnabled(False)

    # Appends to list of moment widgets using default values
    def appendMomentBoxes(self):
        combobox = QComboBox(self.ui.momentTable)
        combobox.addItem("None")
        combobox.addItem("Use seed")
        combobox.addItem("Specify")
        spinbox = QDoubleSpinBox(self.ui.momentTable)
        spinbox.setDecimals(4)
        spinbox.setMaximum(10.0)
        spinbox.setMinimum(-10.0)
        spinbox.setSingleStep(0.01)
        spinbox.setValue(-0.1)
        spinbox.setEnabled(False)
        self.momenttype_boxes.append(combobox)
        self.minmoment_boxes.append(spinbox)
        nmoment = len(self.momenttype_boxes)
        self.momenttype_boxes[nmoment-1].currentIndexChanged[str].connect(
                                         self.setMomentBoxesEnabled)

    # Sets up moment constraints table using default values
    def setupMomentTable(self):
        self.momenttype_boxes = []
        self.minmoment_boxes = []
        for i in range(operatingpoints.numPoints()):
            self.appendMomentBoxes()
            self.ui.momentTable.setRowCount(i+1)
            self.ui.momentTable.setItem(i, 0, QTableWidgetItem("{:d}".format(i+1)))
            self.ui.momentTable.setCellWidget(i, 1, self.momenttype_boxes[i])
            self.ui.momentTable.setCellWidget(i, 2, self.minmoment_boxes[i])
        self.ui.momentTable.verticalHeader().setVisible(False)

    # Sets min moment boxes enabled or disabled
    def setMomentBoxesEnabled(self):
        for i in range(len(self.momenttype_boxes)):
            if self.momenttype_boxes[i].currentText() == "Specify":
                self.minmoment_boxes[i].setEnabled(True)
            else:
                self.minmoment_boxes[i].setEnabled(False)

    # Appends to list of Cp widgets using default values
    def appendCpBoxes(self):
        combobox = QComboBox(self.ui.cpTable)
        combobox.addItem("None")
        combobox.addItem("Specify")
        spinbox = QDoubleSpinBox(self.ui.cpTable)
        spinbox.setDecimals(4)
        spinbox.setMaximum(0.0)
        spinbox.setMinimum(-99.99)
        spinbox.setSingleStep(0.01)
        spinbox.setValue(-2.0)
        spinbox.setEnabled(False)
        self.cptype_boxes.append(combobox)
        self.mincp_boxes.append(spinbox)
        ncp = len(self.cptype_boxes)
        self.cptype_boxes[ncp-1].currentIndexChanged[str].connect(self.setCpBoxesEnabled)

    # Sets up Cp constraints table with default values
    def setupCpTable(self):
        self.cptype_boxes = []
        self.mincp_boxes = []
        for i in range(operatingpoints.numPoints()):
            self.appendCpBoxes()
            self.ui.cpTable.setRowCount(i+1)
            self.ui.cpTable.setItem(i, 0, QTableWidgetItem("{:d}".format(i+1)))
            self.ui.cpTable.setCellWidget(i, 1, self.cptype_boxes[i])
            self.ui.cpTable.setCellWidget(i, 2, self.mincp_boxes[i])
        self.ui.cpTable.verticalHeader().setVisible(False)

    # Sets min Cp boxes enabled or disabled
    def setCpBoxesEnabled(self):
        for i in range(len(self.cptype_boxes)):
            if self.cptype_boxes[i].currentText() == "Specify":
                self.mincp_boxes[i].setEnabled(True)
            else:
                self.mincp_boxes[i].setEnabled(False)

    # Adds a thickness constraint
    def addThicknessConstraint(self):
        xbox = QDoubleSpinBox(self.ui.thicknessTable)
        xbox.setDecimals(4)
        xbox.setMaximum(1.0)
        xbox.setMinimum(0.0)
        xbox.setSingleStep(0.05)
        xbox.setValue(0.5)
        minbox = QDoubleSpinBox(self.ui.thicknessTable)
        minbox.setDecimals(4)
        minbox.setMaximum(99.99)
        minbox.setMinimum(0.0)
        minbox.setSingleStep(1.0)
        minbox.setValue(0.0)
        maxbox = QDoubleSpinBox(self.ui.thicknessTable)
        maxbox.setDecimals(4)
        maxbox.setMaximum(1000.0)
        maxbox.setMinimum(0.0)
        maxbox.setSingleStep(1.0)
        maxbox.setValue(1000.0)
        self.addthickx_boxes.append(xbox)
        self.addthickmin_boxes.append(minbox)
        self.addthickmax_boxes.append(maxbox)

        nrows = self.ui.thicknessTable.rowCount()
        self.ui.thicknessTable.setRowCount(nrows+1)
        self.ui.thicknessTable.setCellWidget(nrows, 0, self.addthickx_boxes[nrows])
        self.ui.thicknessTable.setCellWidget(nrows, 1, self.addthickmin_boxes[nrows])
        self.ui.thicknessTable.setCellWidget(nrows, 2, self.addthickmax_boxes[nrows])
        self.ui.deleteButton.setEnabled(True)

    # Deletes the selected thickness constraint
    def deleteThicknessConstraint(self):
        row = self.ui.thicknessTable.currentRow()
        if row < 0:
            return
        self.ui.thicknessTable.removeRow(row)
        self.addthickx_boxes.pop(row)
        self.addthickmin_boxes.pop(row)
        self.addthickmax_boxes.pop(row)
        if self.ui.thicknessTable.rowCount() == 0:
            self.ui.deleteButton.setEnabled(False)

    def saveConstraints(self):
        constraints.minThickness = self.ui.minThicknessBox.value()/100.
        constraints.maxThickness = self.ui.maxThicknessBox.value()/100.
        constraints.minCamber = self.ui.minCamberBox.value()/100.
        constraints.maxCamber = self.ui.maxCamberBox.value()/100.
        constraints.minTEAngle = self.ui.minTEAngleBox.value()
        constraints.checkCurvature = self.ui.checkCurvatureBox.isChecked()
        constraints.maxReverseTop = self.ui.maxReverseTopBox.value()
        constraints.maxReverseBot = self.ui.maxReverseBottomBox.value()
        constraints.curveThreshold = self.ui.curveThresholdBox.value()
        constraints.symmetrical = self.ui.symmetricalBox.isChecked()
        constraints.minFlapAngle = self.ui.minFlapAngleBox.value()
        constraints.maxFlapAngle = self.ui.maxFlapAngleBox.value()

        # Moment constraints
        constraints.clearMomentConstraints()
        for i in range(len(self.momenttype_boxes)):
            constraint_type = self.momenttype_boxes[i].currentText()
            min_moment = self.minmoment_boxes[i].value()
            constraints.addMomentConstraint(constraint_type, min_moment)

        # Cp constraints
        constraints.clearCpConstraints()
        for i in range(len(self.cptype_boxes)):
            constraint_type = self.cptype_boxes[i].currentText()
            min_cp = self.mincp_boxes[i].value()
            constraints.addCpConstraint(constraint_type, min_cp)

        # Additional thickness constraints
        constraints.clearThicknessConstraints()
        for i in range(len(self.addthickx_boxes)):
            x = self.addthickx_boxes[i].value()
            thickmin = self.addthickmin_boxes[i].value()/100.
            thickmax = self.addthickmax_boxes[i].value()/100
            constraints.addThicknessConstraint(x, thickmin, thickmax)
