#include "settingswindow.h"
#include "settingsbrowser.h"

#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>
#include <QString>
#include <QFont>
#include <QIcon>

/******************************************************************************/
//
// Constructor
//
/******************************************************************************/
SettingsBrowser::SettingsBrowser(SettingsWindow *parent) : QListWidget(parent)
{
  QFont myfont;

  // Store pointer to settings window

  settingswindow = parent;

  // List items
  
  myfont.setBold(true);
  QListWidgetItem *lbl = new QListWidgetItem("Settings Browser", this);
  lbl->setFont(myfont);
  lbl->setIcon(QIcon(":/icons/settings.png"));

  optimitem = new QListWidgetItem("Optimization", this);
  operitem = new QListWidgetItem("Operating conditions", this);
  constritem = new QListWidgetItem("Constraints", this);
  inititem = new QListWidgetItem("Initialization", this);
  psoitem = new QListWidgetItem("Particle swarm", this);
  gaitem = new QListWidgetItem("Genetic algorithm", this);
  simplexitem = new QListWidgetItem("Simplex search", this);
  xfanaitem = new QListWidgetItem("Xfoil analysis", this);
  xfpanitem = new QListWidgetItem("Xfoil paneling", this);

  // Add list items to listwidget

  insertItem(0, lbl);
  insertItem(1, optimitem);
  insertItem(2, operitem);
  insertItem(3, constritem);
  insertItem(4, inititem);
  insertItem(5, psoitem);
  insertItem(6, gaitem);
  insertItem(7, simplexitem);
  insertItem(8, xfanaitem);
  insertItem(9, xfpanitem);

  // Set item properties

  lbl->setFlags(lbl->flags() & ~Qt::ItemIsSelectable);
  optimitem->setCheckState(Qt::Unchecked);
  optimitem->setFlags(optimitem->flags() & ~Qt::ItemIsUserCheckable);
  operitem->setCheckState(Qt::Unchecked);
  operitem->setFlags(operitem->flags() & ~Qt::ItemIsUserCheckable);
  constritem->setCheckState(Qt::Unchecked);
  constritem->setFlags(constritem->flags() & ~Qt::ItemIsUserCheckable);
  inititem->setCheckState(Qt::Checked);
  inititem->setFlags(inititem->flags() & ~Qt::ItemIsUserCheckable);
  psoitem->setCheckState(Qt::Checked);
  psoitem->setFlags(psoitem->flags() & ~Qt::ItemIsUserCheckable);
  gaitem->setCheckState(Qt::Checked);
  gaitem->setFlags(gaitem->flags() & ~Qt::ItemIsUserCheckable);
  simplexitem->setCheckState(Qt::Checked);
  simplexitem->setFlags(simplexitem->flags() & ~Qt::ItemIsUserCheckable);
  xfanaitem->setCheckState(Qt::Checked);
  xfanaitem->setFlags(xfanaitem->flags() & ~Qt::ItemIsUserCheckable);
  xfpanitem->setCheckState(Qt::Checked);
  xfpanitem->setFlags(xfpanitem->flags() & ~Qt::ItemIsUserCheckable);

  // Set tooltips

  optimitem->setToolTip("There are required items left to complete for " +
                        QString("Optimization."));
  operitem->setToolTip("There are required items left to complete for " +
                       QString("Operating conditions."));
  constritem->setToolTip("There are required items left to complete for " +
                         QString("Constraints."));
  inititem->setToolTip("No required items left to complete for " +
                       QString("Initialization."));
  psoitem->setToolTip("No required items left to complete for " +
                      QString("Particle swarm."));
  gaitem->setToolTip("No required items left to complete for " +
                     QString("Genetic algorithm."));
  simplexitem->setToolTip("No required items left to complete for " +
                          QString("Simplex search."));
  xfanaitem->setToolTip("No required items left to complete for " +
                        QString("Xfoil analysis."));
  xfpanitem->setToolTip("No required items left to complete for " +
                        QString("Xfoil paneling."));

  // Connect selection signal/slot

  connect(this, &QListWidget::currentItemChanged, this, 
          &SettingsBrowser::selectionChanged); 

  // Set current item

  setCurrentItem(optimitem);
}

/******************************************************************************/
//
// Signals parent SettingsWindow to change settings pane
//
/******************************************************************************/
void SettingsBrowser::selectionChanged ()
{
  if (currentItem()->text() == "Optimization")
  {
    settingswindow->showOptSettings();
  }
  else if (currentItem()->text() == "Operating conditions")
  {
    settingswindow->showOperSettings();
  }
  else if (currentItem()->text() == "Constraints")
  {
    settingswindow->showConstrSettings();
  }
}
