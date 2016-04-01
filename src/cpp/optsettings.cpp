#include "optsettings.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QSize>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for optimization settings class
//
/******************************************************************************/
OptSettings::OptSettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl, *caselbl, *searchlbl, *globallbl, *locallbl;
  QGridLayout *grid;
  QVBoxLayout *vbox;

  // Label for optimization settings

  lbl = new QLabel("<b>Optimization settings</b>", this);

  // Case name

  caselbl = new QLabel("Case name", this);
  caseedit = new QLineEdit("optfoil", this);

  // Search type

  searchlbl = new QLabel("Search type", this);
  searchbox = new QComboBox(this);
  searchbox->addItem("Global + local");
  searchbox->addItem("Global");
  searchbox->addItem("Local");

  // Global search

  globallbl = new QLabel("Global search", this);
  globalbox = new QComboBox(this);
  globalbox->addItem("Particle swarm");
  globalbox->addItem("Genetic algorithm");
  globalbtn = new QPushButton(QIcon(":/icons/settings.png"), "", this);
  globalbtn->setToolTip("Particle swarm settings");
  globalbtn->setFlat(true);
  globalbtn->setIconSize(QSize(20,20));

  // Local search

  locallbl = new QLabel("Local search", this);
  localbox = new QComboBox(this);
  localbox->addItem("Simplex search");
  localbtn = new QPushButton(QIcon(":/icons/settings.png"), "", this);
  localbtn->setToolTip("Simplex search settings");
  localbtn->setFlat(true);
  localbtn->setIconSize(QSize(20,20));

  // Grid layout
  
  grid = new QGridLayout();
  grid->addWidget(caselbl, 0, 0);
  grid->addWidget(caseedit, 0, 1);
  grid->addWidget(searchlbl, 1, 0);
  grid->addWidget(searchbox, 1, 1);
  grid->addWidget(globallbl, 2, 0);
  grid->addWidget(globalbox, 2, 1);
  grid->addWidget(globalbtn, 2, 2);
  grid->addWidget(locallbl, 3, 0);
  grid->addWidget(localbox, 3, 1);
  grid->addWidget(localbtn, 3, 2);

  // Box layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addLayout(grid);
  vbox->addStretch(0);

  setLayout(vbox);
}
