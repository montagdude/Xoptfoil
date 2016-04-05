#include "settingswindow.h"
#include "optsettings.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QSize>
#include <QString>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for optimization settings class
//
/******************************************************************************/
OptSettings::OptSettings(SettingsWindow *parent) : QWidget(parent)
{
  QGridLayout *grid1, *grid2;
  QHBoxLayout *hbox1, *hbox2, *hbox3, *hbox4;
  QVBoxLayout *vbox;

  // Store pointer to settings window

  settingswindow = parent;

  // Case name

  caseedit = new QLineEdit("optfoil", this);

  // Search type

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

  seedbox = new QComboBox(this);
  seedbox->addItem("From file");
  seedbox->addItem("NACA 4-digit");

  // Select seed airfoil from file

  seedfilelbl = new QLabel("File", this);
  seedfilebox = new QLineEdit(this);
  seedfilebtn = new QPushButton(QIcon(":/icons/browse.png"), " Browse", this);
  seedfilebtn->setFlat(true);
  seedfilebtn->setIconSize(QSize(20,20));

  // NACA 4 digits

  digitlbl = new QLabel("NACA digits", this);
  digitlbl->setEnabled(false);
  digitbox1 = new QSpinBox(this);
  digitbox1->setMinimum(0);
  digitbox1->setValue(0);
  digitbox2 = new QSpinBox(this);
  digitbox2->setMinimum(0);
  digitbox2->setMaximum(9);
  digitbox2->setValue(0);
  digitbox3 = new QSpinBox(this);
  digitbox3->setMinimum(0);
  digitbox3->setValue(1);
  digitbox4 = new QSpinBox(this);
  digitbox4->setMinimum(0);
  digitbox4->setValue(2);

  // Top grid layout (basic settings)
  
  grid1 = new QGridLayout();
  grid1->addWidget(new QLabel("Case name", this), 0, 0);
  grid1->addWidget(caseedit, 0, 1);
  grid1->addWidget(new QLabel("Search type", this), 1, 0);
  grid1->addWidget(searchbox, 1, 1);
  grid1->addWidget(globallbl, 2, 0);
  grid1->addWidget(globalbox, 2, 1);
  hbox1 = new QHBoxLayout();
  hbox1->addWidget(globalbtn);
  hbox1->addStretch(0);
  grid1->addLayout(hbox1, 2, 2);
  grid1->addWidget(locallbl, 3, 0);
  grid1->addWidget(localbox, 3, 1);
  hbox2 = new QHBoxLayout();
  hbox2->addWidget(localbtn);
  hbox2->addStretch(0);
  grid1->addLayout(hbox2, 3, 2);
  grid1->addWidget(new QLabel("Seed airfoil", this), 4, 0);
  grid1->addWidget(seedbox, 4, 1);
  hbox3 = new QHBoxLayout();
  hbox3->addWidget(seedfilelbl);
  hbox3->addWidget(seedfilebox);
  hbox3->addWidget(seedfilebtn);
  grid1->addLayout(hbox3, 4, 2);
  hbox4 = new QHBoxLayout();
  hbox4->addWidget(digitlbl);
  hbox4->addWidget(digitbox1);
  hbox4->addWidget(digitbox2);
  hbox4->addWidget(digitbox3);
  hbox4->addWidget(digitbox4);
  grid1->addLayout(hbox4, 5, 2);

  // Shape functions

  shapebox = new QComboBox(this);
  shapebox->addItem("Hicks-Henne");
  shapebox->addItem("NACA");

  // Number of shape functions

  nshapetbox = new QSpinBox(this);
  nshapebbox = new QSpinBox(this);
  nshapetbox->setMinimum(1);
  nshapetbox->setValue(4);
  nshapebbox->setMinimum(1);
  nshapebbox->setValue(4);

  // Initial perturbation

  initperturbbox = new QDoubleSpinBox(this);
  initperturbbox->setDecimals(3);
  initperturbbox->setValue(0.025);
  initperturbbox->setMinimum(0.0);
  initperturbbox->setSingleStep(0.005);
  initperturbbox->setToolTip("Max change in airfoil surface " +
                             QString("during initialization"));

  // Minimum bump width

  minbumpbox = new QDoubleSpinBox(this);
  minbumpbox->setValue(0.1);
  minbumpbox->setSingleStep(0.1);
  minbumpbox->setMinimum(0.0);
  minbumpbox->setMaximum(1.0);
  minbumpbox->setDecimals(2);

  // Restart write frequency

  restbox = new QSpinBox(this);
  restbox->setMinimum(1);
  restbox->setValue(20);

  // Whether to write design progression

  writedesignbox = new QCheckBox("Write design progression to files", this);
  writedesignbox->setChecked(true);

  // Bottom grid layout (advanced settings)

  grid2 = new QGridLayout();
  grid2->addWidget(new QLabel("Shape functions", this), 0, 0);
  grid2->addWidget(shapebox, 0, 1);
  grid2->addWidget(new QLabel("Number of functions for top surface", this), 
                   1, 0);
  grid2->addWidget(nshapetbox, 1, 1);
  grid2->addWidget(new QLabel("Number of functions for bottom surface", this),
                   2, 0);
  grid2->addWidget(nshapebbox, 2, 1);
  grid2->addWidget(new QLabel("Initial perturbation", this), 3, 0);
  grid2->addWidget(initperturbbox, 3, 1);
  grid2->addWidget(new QLabel("Minimum bump width", this), 4, 0);
  grid2->addWidget(minbumpbox, 4, 1);
  grid2->addWidget(new QLabel("Restart write frequency", this), 5, 0);
  grid2->addWidget(restbox, 5, 1);
  grid2->addWidget(writedesignbox, 6, 0);

  // Box layout

  vbox = new QVBoxLayout();
  vbox->addWidget(new QLabel("<b>Basic optimization settings</b>", this));
  vbox->addLayout(grid1);
  vbox->addWidget(new QLabel("<b>Advanced settings</b>", this));
  vbox->addLayout(grid2);
  vbox->addStretch(0);

  setLayout(vbox);

  // Connect signals and slots (note: SIGNAL/SLOT syntax seems necessary for
  // signals with an argument; also needs Q_OBJECT in header)

  connect(searchbox, SIGNAL(currentIndexChanged(int)), this,
          SLOT(searchBoxChanged(int)));
  connect(globalbox, SIGNAL(currentIndexChanged(int)), this,
          SLOT(globalBoxChanged(int)));
  connect(globalbtn, &QPushButton::clicked, this, 
          &OptSettings::globalBtnClicked);
  connect(localbtn, &QPushButton::clicked, this, 
          &OptSettings::localBtnClicked);
  connect(seedbox, SIGNAL(currentIndexChanged(int)), this, 
          SLOT(seedBoxChanged(int)));
}

/******************************************************************************/
//
// Defines behavior when search type box index is changed
// 
/******************************************************************************/
void OptSettings::searchBoxChanged ( int idx )
{
  // Global + local
  
  if (idx == 0)
  {
    globallbl->setEnabled(true);
    globalbox->setEnabled(true);
    globalbtn->setEnabled(true);
    locallbl->setEnabled(true);
    localbox->setEnabled(true);
    localbtn->setEnabled(true);
  }

  // Global
  
  else if (idx == 1)
  {
    globallbl->setEnabled(true);
    globalbox->setEnabled(true);
    globalbtn->setEnabled(true);
    locallbl->setEnabled(false);
    localbox->setEnabled(false);
    localbtn->setEnabled(false);
  }

  // Local
  
  else if (idx == 2)
  {
    globallbl->setEnabled(false);
    globalbox->setEnabled(false);
    globalbtn->setEnabled(false);
    locallbl->setEnabled(true);
    localbox->setEnabled(true);
    localbtn->setEnabled(true);
  }
}

/******************************************************************************/
//
// Defines behavior when global search box index is changed
//
/******************************************************************************/
void OptSettings::globalBoxChanged ( int idx )
{
  // Particle swarm

  if (idx == 0) { globalbtn->setText("Particle swarm settings"); }

  // Genetic algorithm
  
  else if (idx == 1) { globalbtn->setText("Genetic algorithm settings"); }
}

/******************************************************************************/
//
// Defines behavior when global search settings button is clicked
//
/******************************************************************************/
void OptSettings::globalBtnClicked ()
{
  // Particle swarm settings

  if (globalbox->currentIndex() == 0) { settingswindow->showPSOSettings(); }

  // Genetic algorithm settings

  else if (globalbox->currentIndex() == 1) { settingswindow->showGASettings(); }
}

/******************************************************************************/
//
// Defines behavior when local search settings button is clicked
//
/******************************************************************************/
void OptSettings::localBtnClicked ()
{
  // Simplex search settings

  if (localbox->currentIndex() == 0) { settingswindow->showSimplexSettings();}
}

/******************************************************************************/
//
// Defines behavior when seed airfoil selection box index is changed
//
/******************************************************************************/
void OptSettings::seedBoxChanged ( int idx )
{
  // Seed airfoil from file
  
  if (idx == 0)
  {
    digitlbl->setEnabled(false);
    digitbox1->setEnabled(false);
    digitbox2->setEnabled(false);
    digitbox3->setEnabled(false);
    digitbox4->setEnabled(false);
    seedfilelbl->setEnabled(true);
    seedfilebox->setEnabled(true);
    seedfilebtn->setEnabled(true);
  }

  // NACA 4-digit
  
  else if (idx == 1)
  {
    seedfilelbl->setEnabled(false);
    seedfilebox->setEnabled(false);
    seedfilebtn->setEnabled(false);
    digitlbl->setEnabled(true);
    digitbox1->setEnabled(true);
    digitbox2->setEnabled(true);
    digitbox3->setEnabled(true);
    digitbox4->setEnabled(true);
  }
}
