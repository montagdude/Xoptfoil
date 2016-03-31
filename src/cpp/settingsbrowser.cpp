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
SettingsBrowser::SettingsBrowser(QWidget *parent) : QListWidget(parent)
{
  QFont myfont;

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
  xfrunitem = new QListWidgetItem("Xfoil analysis", this);
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
  insertItem(8, xfrunitem);
  insertItem(9, xfpanitem);

  // Set item properties

  lbl->setFlags(Qt::ItemIsEnabled);
  optimitem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled); 
  optimitem->setCheckState(Qt::Unchecked);
  operitem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled); 
  operitem->setCheckState(Qt::Unchecked);
  constritem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled); 
  constritem->setCheckState(Qt::Unchecked);

  // Set tooltips

  optimitem->setToolTip("There are items left to complete for Optimization.");
  operitem->setToolTip("There are items left to complete for Operating " +
                       QString("conditions."));
  constritem->setToolTip("There are items left to complete for Constraints.");

  // Set current item

  setCurrentItem(optimitem);
}
