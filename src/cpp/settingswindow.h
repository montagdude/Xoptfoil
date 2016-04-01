#pragma once

#include <QWidget>
#include <QSplitter>
#include <QStackedWidget>
#include <QString>

// Forward declarations

class SettingsBrowser;
class OptSettings;
class OperSettings;

/******************************************************************************/
//
// Header for settings window class
//
/******************************************************************************/
class SettingsWindow : public QSplitter
{
  private:

    // Widgets contained in splitter

    OptSettings *opt_settings;
    OperSettings *oper_settings;
    QStackedWidget *settings_pane;
    SettingsBrowser *settingsbrowser;

  public:

    // Constructor

    SettingsWindow ( QWidget *parent = 0 );

  public slots:

    // Showing different settings panes

    void showOptSettings ();
    void showOperSettings ();
};
