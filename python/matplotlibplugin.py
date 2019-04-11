import os
from PyQt5.QtGui import QIcon
from PyQt5.QtDesigner import QPyDesignerCustomWidgetPlugin
from matplotlib import rcParams
from matplotlibwidget import MatplotlibWidget

rcParams['font.size'] = 9


class MatplotlibPlugin(QPyDesignerCustomWidgetPlugin):
    def __init__(self, parent=None):
        super(MatplotlibPlugin, self).__init__(parent)
        self._initialized = False

    def initialize(self, editor):
        self._initialized = True

    def isInitialized(self):
        return self._initialized

    def createWidget(self, parent):
        return MatplotlibWidget(parent)

    def name(self):
        return 'MatplotlibWidget'

    def group(self):
        return 'PyQt'

    def icon(self):
        return QIcon(os.path.join(
            rcParams['datapath'], 'images', 'matplotlib.png'))

    def toolTip(self):
        return ''

    def whatsThis(self):
        return ''

    def isContainer(self):
        return False

    def domXml(self):
        return '<widget class="MatplotlibWidget" name="mplwidget">\n' \
               '</widget>\n'

    def includeFile(self):
        return 'matplotlibwidget'
