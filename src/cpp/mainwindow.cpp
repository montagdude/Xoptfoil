#include "settingsbrowser.h"
#include "mainwindow.h"

#include <QWidget>
#include <QSplitter>
#include <QScrollArea>
#include <QList>
#include <QMainWindow>

/******************************************************************************/
//
// Constructor for main window
//
/******************************************************************************/
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
  QScrollArea *scrollarea;
  QSplitter *splitter;
  QList<int> default_sizes;

  // Settings browser

  settingsbrowser = new SettingsBrowser(this);

  // Scroll area

  scrollarea = new QScrollArea(this);

  // Splitter (allows user to resize widgets)

  splitter = new QSplitter(this);
  splitter->setOrientation(Qt::Horizontal);
  splitter->addWidget(settingsbrowser);
  splitter->addWidget(scrollarea);
  
  // Default sizes

  default_sizes.append(220);
  default_sizes.append(680); 
  splitter->setSizes(default_sizes);

  // Central widget
  
  setCentralWidget(splitter);
}
