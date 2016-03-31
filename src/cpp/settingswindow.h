#pragma once

#include "settingsbrowser.h"

#include <QWidget>
#include <QSplitter>
#include <QScrollArea>

/******************************************************************************/
//
// Header for settings window class
//
/******************************************************************************/
class SettingsWindow : public QSplitter
{
  private:

    // Widgets contained in splitter

    SettingsBrowser *settingsbrowser;
    QScrollArea *scrollarea;
    
  public:

    // Constructor

    SettingsWindow ( QWidget *parent = 0 );
};
