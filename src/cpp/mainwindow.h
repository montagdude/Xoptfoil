#pragma once

#include "navigatorlist.h"

#include <QWidget>
#include <QMainWindow>

/******************************************************************************/
//
// Header for main window
//
/******************************************************************************/
class MainWindow : public QMainWindow
{
  private:

    // Items in main window

    NavigatorList *navigatorlist;

  public:
 
    // Constructor

    MainWindow ( QWidget *parent = 0 );
};
