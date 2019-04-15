from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QWidget, QSizePolicy
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlibwidget import MatplotlibWidget

class MatplotlibToolbar(NavigationToolbar):
    def __init__(self, parent=None, canvas=None):
        # If canvas is not supplied, set a temporary one. In that case, you should connect it
        # later with setCanvas.
        if canvas is None:
            canvas = MatplotlibWidget()

        super(MatplotlibToolbar, self).__init__(canvas, QWidget())
        self.setParent(parent)
        super(MatplotlibToolbar, self).setSizePolicy(QSizePolicy.Expanding,
                                                     QSizePolicy.Minimum)

    def minimumSizeHint(self):
        return QSize(10,5)

    def setCanvas(self, canvas):
        self.canvas = canvas
