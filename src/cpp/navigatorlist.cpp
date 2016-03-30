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

  globalitem = new QListWidgetItem("Optimization settings", this);
  operitem = new QListWidgetItem("Operating conditions", this);
  constritem = new QListWidgetItem("Constraints", this);
  inititem = new QListWidgetItem("Initialization", this);
  psoitem = new QListWidgetItem("Particle swarm settings", this);
  gaitem = new QListWidgetItem("Genetic algorithm settings", this);
  simplexitem = new QListWidgetItem("Simplex search settings", this);
  xfrunitem = new QListWidgetItem("Xfoil run settings", this);
  xfpanitem = new QListWidgetItem("Xfoil paneling settings", this);
}
