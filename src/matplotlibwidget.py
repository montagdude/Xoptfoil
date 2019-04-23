from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure
import numpy as np
from math import pi

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

        #FIXME make these configurable options stored in a startup file
        self.bgcolor = '#000000' 
        self.fgcolor = '#ffffff'
        self.seedcolor = '#4b97ff'
        self.currentcolor = '#ff2424'
        self.figure.set_facecolor(self.bgcolor)
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
        self.ax.clear()
        self.ax.set_facecolor(self.bgcolor)
        self.ax.spines['top'].set_color(self.fgcolor)
        self.ax.spines['bottom'].set_color(self.fgcolor)
        self.ax.spines['left'].set_color(self.fgcolor)
        self.ax.spines['right'].set_color(self.fgcolor)
        self.ax.xaxis.label.set_color(self.fgcolor)
        self.ax.yaxis.label.set_color(self.fgcolor)
        self.ax.tick_params(axis='x', colors=self.fgcolor)
        self.ax.tick_params(axis='y', colors=self.fgcolor)

    def plotAirfoils(self, seed, newdesign=None):
        """Plots seed airfoil and, optionally, new design

        Inputs:
            seed: seed airfoil
            newdesign: current design
        """
        self.setupAxes()
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.plot(seed.x, seed.y, color=self.seedcolor, label='Seed airfoil')
        if newdesign is not None:
            self.ax.plot(newdesign.x, newdesign.y, color=self.currentcolor,
                         label='Current design')
        self.ax.set_aspect('equal', 'datalim')
        leg = self.ax.legend(facecolor=self.bgcolor, edgecolor=self.fgcolor)
        for text in leg.get_texts():
            text.set_color(self.fgcolor)
        self.draw()
