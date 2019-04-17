#!/usr/bin/python3

import sys
import os
from PyQt5 import QtWidgets

# FIXME: this will have to be set at install time. Also doesn't work if not
# run in the src directory, so obviously this is only temporary.
installdir = os.path.join(os.getcwd(), '..')
sys.path.append(os.path.join(installdir, 'ui'))

import mainwindow
from xfoilpaneling_dialog import XfoilPanelingDialog
from data import Data

class XoptfoilMainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(XoptfoilMainWindow, self).__init__()
        self.data = Data()
        self.currpath = os.getcwd()

        self.ui = mainwindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle("Xoptfoil")

        # Connect toolbar to canvas
        self.ui.mpltoolbar.setCanvas(self.ui.mplwidget)

        # Signals and slots
        self.ui.action_Load_seed_airfoil.triggered.connect(self.loadSeed)
        self.ui.action_Quit.triggered.connect(self.close)
        self.ui.action_Xfoil_paneling.triggered.connect(self.showXfoilPaneling)

    # Loads seed airfoil from file
    def loadSeed(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load seed airfoil",
                                                         self.currpath)
        if fname == "":
            return
        else:
            self.currpath = os.path.dirname(fname)

        retval = self.data.readSeedAirfoil(fname)
        errmsg = None
        if retval == 1:
            errmsg  = "Unable to read airfoil file {:s}: I/O error.".format(fname)
        elif retval == 2:
            errmsg  = "Unable to read airfoil file {:s}:".format(fname)
            errmsg += " Not a valid airfoil file."
        if errmsg is not None:
            QtWidgets.QMessageBox.critical(self, "Error", errmsg)
        self.ui.mplwidget.plotAirfoils(self.data.seed_airfoil)

    def showXfoilPaneling(self):
        dialog = XfoilPanelingDialog()
        retval = dialog.exec()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    application = XoptfoilMainWindow()
    application.show()
    sys.exit(app.exec())
