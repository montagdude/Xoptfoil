from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure
import numpy as np
from math import pi

from settings import plotsettings
from data import data

class MatplotlibWidget(Canvas):
    def __init__(self, parent=None, title='', xlabel='', ylabel='',
                 xlim=None, ylim=None, xscale='linear', yscale='linear',
                 width=4, height=3, dpi=100):
        self.figure = Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title(title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        if xscale is not None:
            self.ax.set_xscale(xscale)
        if yscale is not None:
            self.ax.set_yscale(yscale)
        if xlim is not None:
            self.ax.set_xlim(*xlim)
        if ylim is not None:
            self.ax.set_ylim(*ylim)
        self.setupAxes()

        super(MatplotlibWidget, self).__init__(self.figure)
        self.setParent(parent)
        super(MatplotlibWidget, self).setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Expanding)
        super(MatplotlibWidget, self).updateGeometry()

    def sizeHint(self):
        return QSize(*self.get_width_height())

    def minimumSizeHint(self):
        return QSize(10, 10)

    def setupAxes(self):
        self.figure.set_facecolor(plotsettings.value("bgColor"))
        self.ax.clear()
        self.ax.set_facecolor(plotsettings.value("bgColor"))
        self.ax.spines['top'].set_color(plotsettings.value("fgColor"))
        self.ax.spines['bottom'].set_color(plotsettings.value("fgColor"))
        self.ax.spines['left'].set_color(plotsettings.value("fgColor"))
        self.ax.spines['right'].set_color(plotsettings.value("fgColor"))
        self.ax.xaxis.label.set_color(plotsettings.value("fgColor"))
        self.ax.yaxis.label.set_color(plotsettings.value("fgColor"))
        self.ax.tick_params(axis='x', colors=plotsettings.value("fgColor"))
        self.ax.tick_params(axis='y', colors=plotsettings.value("fgColor"))

    def plotAirfoils(self):
        """Plots seed airfoil and, optionally, new design
        """
        self.setupAxes()
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.plot(data.seed_airfoil.x, data.seed_airfoil.y, 
                     color=plotsettings.value("seedColor"), label='Seed airfoil')
        if data.current_airfoil is not None:
            self.ax.plot(data.current_airfoil.x, data.current_airfoil.y,
                         color=plotsettings.value("currentColor"), label='Current design')
        self.ax.set_aspect('equal', 'datalim')
        leg = self.ax.legend(facecolor=plotsettings.value("bgColor"),
                             edgecolor=plotsettings.value("fgColor"))
        for text in leg.get_texts():
            text.set_color(plotsettings.value("fgColor"))
        self.draw()
