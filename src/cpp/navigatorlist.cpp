#include "navigatorlist.h"

#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>

/******************************************************************************/
//
// Constructor
//
/******************************************************************************/
NavigatorList::NavigatorList(QWidget *parent) : QListWidget(parent)
{
  // List items
  
  QListWidgetItem *lbl = new QListWidgetItem("Settings Browser", this);

  optimitem = new QListWidgetItem("Optimization settings", this);
  operitem = new QListWidgetItem("Operating conditions", this);
  constritem = new QListWidgetItem("Constraints", this);
  inititem = new QListWidgetItem("Initialization", this);
  psoitem = new QListWidgetItem("Particle swarm settings", this);
  gaitem = new QListWidgetItem("Genetic algorithm settings", this);
  simplexitem = new QListWidgetItem("Simplex search settings", this);
  xfrunitem = new QListWidgetItem("Xfoil run settings", this);
  xfpanitem = new QListWidgetItem("Xfoil paneling settings", this);

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

  // Set current item

  setCurrentItem(optimitem);
}
