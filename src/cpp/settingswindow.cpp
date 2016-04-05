#include "optsettings.h"
#include "opersettings.h"
#include "constrsettings.h"
#include "initsettings.h"
#include "psosettings.h"
#include "gasettings.h"
#include "simplexsettings.h"
#include "xfanasettings.h"
#include "xfpansettings.h"
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
  init_settings = new QScrollArea(this);
  init_settings->setWidget(new InitSettings(this));
  pso_settings = new QScrollArea(this);
  pso_settings->setWidget(new PSOSettings(this));
  ga_settings = new QScrollArea(this);
  ga_settings->setWidget(new GASettings(this));
  simplex_settings = new QScrollArea(this);
  simplex_settings->setWidget(new SimplexSettings(this));
  xfana_settings = new QScrollArea(this);
  xfana_settings->setWidget(new XfAnaSettings(this));
  xfpan_settings = new QScrollArea(this);
  xfpan_settings->setWidget(new XfPanSettings(this));
  settings_pane = new QStackedWidget(this);
  settings_pane->addWidget(opt_settings);
  settings_pane->addWidget(oper_settings);
  settings_pane->addWidget(constr_settings);
  settings_pane->addWidget(init_settings);
  settings_pane->addWidget(pso_settings);
  settings_pane->addWidget(ga_settings);
  settings_pane->addWidget(simplex_settings);
  settings_pane->addWidget(xfana_settings);
  settings_pane->addWidget(xfpan_settings);
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
  settingsbrowser->setCurrentRow(1);
}

/******************************************************************************/
//
// Switches to operating conditions settings pane
//
/******************************************************************************/
void SettingsWindow::showOperSettings ()
{
  settings_pane->setCurrentWidget(oper_settings);
  settingsbrowser->setCurrentRow(2);
}

/******************************************************************************/
//
// Switches to constraints settings pane
//
/******************************************************************************/
void SettingsWindow::showConstrSettings ()
{
  settings_pane->setCurrentWidget(constr_settings);
  settingsbrowser->setCurrentRow(3);
}

/******************************************************************************/
//
// Switches to initialization settings pane
//
/******************************************************************************/
void SettingsWindow::showInitSettings ()
{
  settings_pane->setCurrentWidget(init_settings);
  settingsbrowser->setCurrentRow(4);
}

/******************************************************************************/
//
// Switches to particle swarm settings pane
//
/******************************************************************************/
void SettingsWindow::showPSOSettings ()
{
  settings_pane->setCurrentWidget(pso_settings);
  settingsbrowser->setCurrentRow(5);
}

/******************************************************************************/
//
// Switches to genetic algorithm settings pane
//
/******************************************************************************/
void SettingsWindow::showGASettings ()
{
  settings_pane->setCurrentWidget(ga_settings);
  settingsbrowser->setCurrentRow(6);
}

/******************************************************************************/
//
// Switches to simplex settings pane
//
/******************************************************************************/
void SettingsWindow::showSimplexSettings ()
{
  settings_pane->setCurrentWidget(simplex_settings);
  settingsbrowser->setCurrentRow(7);
}

/******************************************************************************/
//
// Switches to Xfoil analysis settings pane
//
/******************************************************************************/
void SettingsWindow::showXfAnaSettings ()
{
  settings_pane->setCurrentWidget(xfana_settings);
  settingsbrowser->setCurrentRow(8);
}

/******************************************************************************/
//
// Switches to Xfoil paneling settings pane
//
/******************************************************************************/
void SettingsWindow::showXfPanSettings ()
{
  settings_pane->setCurrentWidget(xfpan_settings);
  settingsbrowser->setCurrentRow(9);
}
