#pragma once

#include "settingsbrowser.h"

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

    SettingsBrowser *settingsbrowser;

  public:
 
    // Constructor

    MainWindow ( QWidget *parent = 0 );
};
