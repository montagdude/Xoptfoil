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
  QLabel *lbl, *caselbl, *searchlbl, *globallbl, *locallbl, *seedlbl;
  QGridLayout *grid;
  QHBoxLayout *hbox1, *hbox2, *hbox3, *hbox4;
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
  globalbtn = new QPushButton(QIcon(":/icons/go.png"), 
                              " Particle swarm settings", this);
  globalbtn->setFlat(true);
  globalbtn->setIconSize(QSize(20,20));

  // Local search

  locallbl = new QLabel("Local search", this);
  localbox = new QComboBox(this);
  localbox->addItem("Simplex search");
  localbtn = new QPushButton(QIcon(":/icons/go.png"), 
                             " Simplex search settings", this);
  localbtn->setFlat(true);
  localbtn->setIconSize(QSize(20,20));

  // Seed airfoil selection

  seedlbl = new QLabel("Seed airfoil", this);
  seedbox = new QComboBox(this);
  seedbox->addItem("From file");
  seedbox->addItem("NACA 4-digit");

  // Select seed airfoil from file

  seedfilebtn = new QPushButton(QIcon(":/icons/browse.png"), " Browse", this);
  seedfilebtn->setFlat(true);
  seedfilebtn->setIconSize(QSize(20,20));

  // NACA 4 digits

  digitlbl = new QLabel("NACA digits", this);
  digitlbl->setEnabled(false);
  digitedit = new QLineEdit("0012");
  digitedit->setMaxLength(4);
  digitedit->setEnabled(false);

  // Grid layout
  
  grid = new QGridLayout();
  grid->addWidget(caselbl, 0, 0);
  grid->addWidget(caseedit, 0, 1);
  grid->addWidget(searchlbl, 1, 0);
  grid->addWidget(searchbox, 1, 1);
  grid->addWidget(globallbl, 2, 0);
  grid->addWidget(globalbox, 2, 1);
  hbox1 = new QHBoxLayout();
  hbox1->addWidget(globalbtn);
  hbox1->addStretch(0);
  grid->addLayout(hbox1, 2, 2);
  grid->addWidget(locallbl, 3, 0);
  grid->addWidget(localbox, 3, 1);
  hbox2 = new QHBoxLayout();
  hbox2->addWidget(localbtn);
  hbox2->addStretch(0);
  grid->addLayout(hbox2, 3, 2);
  grid->addWidget(seedlbl, 4, 0);
  grid->addWidget(seedbox, 4, 1);
  hbox3 = new QHBoxLayout();
  hbox3->addWidget(seedfilebtn);
  hbox3->addStretch(0);
  grid->addLayout(hbox3, 4, 2);
  hbox4 = new QHBoxLayout();
  hbox4->addWidget(digitlbl);
  hbox4->addWidget(digitedit);
  hbox4->addStretch(0);
  grid->addLayout(hbox4, 5, 2);

  // Box layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addLayout(grid);
  vbox->addStretch(0);

  setLayout(vbox);

  // Connect signals and slots (note: this method seems necessary for
  // signals with an argument; also needs Q_OBJECT in header)

  connect(seedbox, SIGNAL(currentIndexChanged(int)), this, 
          SLOT(seedBoxChanged(int)));
}

/******************************************************************************/
//
// Defines behavior when seed airfoil selection box index is changed
//
/******************************************************************************/
void OptSettings::seedBoxChanged ( int idx )
{
  if (idx == 0)
  {
    digitlbl->setEnabled(false);
    digitedit->setEnabled(false);
    seedfilebtn->setEnabled(true);
  }
  else if (idx == 1)
  {
    seedfilebtn->setEnabled(false);
    digitlbl->setEnabled(true);
    digitedit->setEnabled(true);
  }
}
