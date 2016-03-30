#include "navigatorlist.h"
#include "mainwindow.h"

#include <QWidget>
#include <QMainWindow>

/******************************************************************************/
//
// Constructor for main window
//
/******************************************************************************/
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
  // Navigator list

  navigatorlist = new NavigatorList(this);

  // Central widget

  setCentralWidget(navigatorlist);
}
