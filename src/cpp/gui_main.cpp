#include "mainwindow.h"

#include <QApplication>

/******************************************************************************/
//
// Main application
//
/******************************************************************************/
int main(int argc, char *argv[])
{
  QApplication xoptfoil(argc, argv);

  // Create main window

  MainWindow window;
  window.setWindowTitle("XoptFoil");
  window.resize(800,600);
  window.show();

  return xoptfoil.exec();
}
