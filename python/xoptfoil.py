#!/usr/bin/python3

import sys
import os
from PyQt5 import QtWidgets

from mainwindow import Ui_MainWindow
from data import Data

class mywindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(mywindow, self).__init__()
        self.data = Data()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle("Xoptfoil")

        # Signals and slots
        self.ui.action_Load_seed_airfoil.triggered.connect(self.loadSeed)
        self.ui.action_Quit.triggered.connect(self.close)

    # Loads seed airfoil from file
    def loadSeed(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load seed airfoil", os.getcwd())
        if not self.data.readSeedAirfoil(fname):
            msg  = "Unable to read airfoil file {:s}.".format(fname)
            msg += " Not a valid airfoil file."
            QtWidgets.QMessageBox.critical(self, "Error", msg)

if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    application = mywindow()
    application.show()
    sys.exit(app.exec())
