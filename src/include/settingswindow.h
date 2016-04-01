#pragma once

#include <QWidget>
#include <QSplitter>
#include <QStackedWidget>
#include <QScrollArea>
#include <QString>

// Forward declarations

class SettingsBrowser;

/******************************************************************************/
//
// Header for settings window class
//
/******************************************************************************/
class SettingsWindow : public QSplitter
{
  private:

    // Widgets contained in splitter

    QScrollArea *opt_settings;
    QScrollArea *oper_settings;
    QScrollArea *constr_settings;
    QScrollArea *init_settings;
    QScrollArea *pso_settings;
    QScrollArea *ga_settings;
    QScrollArea *simplex_settings;
    QScrollArea *xfana_settings;
    QScrollArea *xfpan_settings;
    QStackedWidget *settings_pane;
    SettingsBrowser *settingsbrowser;

  public:

    // Constructor

    SettingsWindow ( QWidget *parent = 0 );

  public slots:

    // Showing different settings panes

    void showOptSettings ();
    void showOperSettings ();
    void showConstrSettings ();
    void showInitSettings ();
    void showPSOSettings ();
    void showGASettings ();
    void showSimplexSettings ();
    void showXfAnaSettings ();
    void showXfPanSettings ();
};
