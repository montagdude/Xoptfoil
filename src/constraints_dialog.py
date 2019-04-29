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
        self.ui.addButton.clicked.connect(self.addThicknessConstraint)
        self.ui.deleteButton.clicked.connect(self.deleteThicknessConstraint)

        # Sets up moment constraints table with default values
        self.momenttype_boxes = []
        self.minmoment_boxes = []
        self.setupMomentTable()

        # Additional thickness constraints widgets
        self.addthickx_boxes = []
        self.addthickmin_boxes = []
        self.addthickmax_boxes = []

    # Sets up moment constraints table with default values
    def setupMomentTable(self):
        self.momenttype_boxes = []
        self.minmoment_boxes = []
        for i in range(operatingpoints.numPoints()):
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
            self.momenttype_boxes[i].currentIndexChanged[str].connect(
                                     self.setMomentBoxesEnabled)
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
        pass
