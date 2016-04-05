#include "settingsbrowser.h"
#include "settingswindow.h"
#include "mainwindow.h"

#include <QWidget>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QApplication>

/******************************************************************************/
//
// Constructor for main window
//
/******************************************************************************/
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
  QMenu *filemenu, *airfoilsmenu, *settingsmenu, *optimizemenu;
  QAction *workdiract, *saveact, *savenmlact, *openact, *readnmlact, *quitact;
  QAction *gotosettingsact, *optimact, *operact, *constract, *initact, *psoact,
          *gaact, *simplexact, *xfanaact, *xfpanact;
  QAction *managefoilsact;
  QAction *gotooptimact, *beginoptimact, *pauseoptimact, *stopoptimact;

  // Create settings window and set as central widget

  settingswindow = new SettingsWindow(this);
  setCentralWidget(settingswindow);

  // Actions for file menu

  workdiract = new QAction("&Set working directory", this);
  saveact = new QAction("Save &case", this);
  openact = new QAction("&Open case", this);
  savenmlact = new QAction("Save settings to &namelist file", this);
  readnmlact = new QAction("&Read settings from namelist file", this);
  quitact = new QAction("&Quit", this);

  // File menu

  filemenu = menuBar()->addMenu("&File");
  filemenu->addAction(workdiract);
  filemenu->addAction(saveact);
  filemenu->addAction(openact);
  filemenu->addAction(savenmlact);
  filemenu->addAction(readnmlact);
  filemenu->addSeparator();
  filemenu->addAction(quitact);

  // Actions for settings menu

  gotosettingsact = new QAction("&Go to settings window ...", this);
  gotosettingsact->setEnabled(false);
  optimact = new QAction("&Optimization", this);
  operact = new QAction("Op&erating conditions", this);
  constract = new QAction("&Constraints", this);
  initact = new QAction("&Initialization", this);
  psoact = new QAction("&Particle swarm", this);
  gaact = new QAction("&Genetic algorithm", this); 
  simplexact = new QAction("Si&mplex search", this);
  xfanaact = new QAction("&Xfoil analysis", this);
  xfpanact = new QAction("X&foil paneling", this);

  // Settings menu

  settingsmenu = menuBar()->addMenu("&Settings");
  settingsmenu->addAction(gotosettingsact);
  settingsmenu->addAction(optimact);
  settingsmenu->addAction(operact);
  settingsmenu->addAction(constract);
  settingsmenu->addAction(initact);
  settingsmenu->addAction(psoact);
  settingsmenu->addAction(gaact);
  settingsmenu->addAction(simplexact);
  settingsmenu->addAction(xfanaact);
  settingsmenu->addAction(xfpanact);

  // Actions for optimization menu
 
  gotooptimact = new QAction("&Go to optimization window ...", this);
  gotooptimact->setEnabled(false);
  beginoptimact = new QAction("&Begin optimization", this);
  beginoptimact->setEnabled(false);
  pauseoptimact = new QAction("&Pause optimization", this);
  pauseoptimact->setEnabled(false);
  stopoptimact = new QAction("&Stop optimization", this);
  stopoptimact->setEnabled(false);

  // Optimization menu
   
  optimizemenu = menuBar()->addMenu("&Optimization");
  optimizemenu->addAction(gotooptimact);
  optimizemenu->addAction(beginoptimact);
  optimizemenu->addAction(pauseoptimact);
  optimizemenu->addAction(stopoptimact);

  // Actions for airfoils menu

  managefoilsact = new QAction("&Manage seed airfoils ...", this);

  // Airfoils menu

  airfoilsmenu = menuBar()->addMenu("&Airfoils");
  airfoilsmenu->addAction(managefoilsact);

  // Connect signals/slots

  connect(quitact, &QAction::triggered, qApp, &QApplication::quit);
  connect(optimact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showOptSettings);
  connect(operact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showOperSettings);
  connect(constract, &QAction::triggered, settingswindow, 
          &SettingsWindow::showConstrSettings);
  connect(initact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showInitSettings);
  connect(psoact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showPSOSettings);
  connect(gaact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showGASettings);
  connect(simplexact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showSimplexSettings);
  connect(xfanaact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showXfAnaSettings);
  connect(xfpanact, &QAction::triggered, settingswindow, 
          &SettingsWindow::showXfPanSettings);
}
