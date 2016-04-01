#include "optsettings.h"
#include "opersettings.h"
#include "constrsettings.h"
#include "settingsbrowser.h"
#include "settingswindow.h"

#include <QWidget>
#include <QSplitter>
#include <QStackedWidget>
#include <QScrollArea>
#include <QList>
#include <QString>

/******************************************************************************/
//
// Constructor
//
/******************************************************************************/
SettingsWindow::SettingsWindow(QWidget *parent) : QSplitter(parent)
{
  QList<int> default_sizes;

  // Widgets contained in splitter
  
  opt_settings = new QScrollArea(this);
  opt_settings->setWidget(new OptSettings(this));
  oper_settings = new QScrollArea(this);
  oper_settings->setWidget(new OperSettings(this));
  constr_settings = new QScrollArea(this);
  constr_settings->setWidget(new ConstrSettings(this));
  settings_pane = new QStackedWidget(this);
  settings_pane->addWidget(opt_settings);
  settings_pane->addWidget(oper_settings);
  settings_pane->addWidget(constr_settings);
  settingsbrowser = new SettingsBrowser(this);

  // Add widgets to splitter and set orientation

  addWidget(settingsbrowser);
  addWidget(settings_pane);
  setOrientation(Qt::Horizontal);

  // Default sizes for widgets
  
  default_sizes.append(220);
  default_sizes.append(680);
  setSizes(default_sizes);
}

/******************************************************************************/
//
// Switches to optimization settings pane
//
/******************************************************************************/
void SettingsWindow::showOptSettings ()
{
  settings_pane->setCurrentWidget(opt_settings);
}

/******************************************************************************/
//
// Switches to operating conditions settings pane
//
/******************************************************************************/
void SettingsWindow::showOperSettings ()
{
  settings_pane->setCurrentWidget(oper_settings);
}

/******************************************************************************/
//
// Switches to constraints settings pane
//
/******************************************************************************/
void SettingsWindow::showConstrSettings ()
{
  settings_pane->setCurrentWidget(constr_settings);
}
