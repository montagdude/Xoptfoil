#include "optsettings.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QSpinBox>
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
OptSettings::OptSettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl, *caselbl, *searchlbl, *globallbl, *locallbl, *seedlbl;
  QLabel *advlbl, *shapelbl, *nshapetlbl, *nshapeblbl, *initperturblbl;
  QLabel *minbumplbl, *restlbl;
  QGridLayout *grid1, *grid2;
  QHBoxLayout *hbox1, *hbox2, *hbox3, *hbox4;
  QVBoxLayout *vbox;

  // Label for optimization settings

  lbl = new QLabel("<b>Basic optimization settings</b>", this);

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

  seedfilelbl = new QLabel("File", this);
  seedfilebox = new QLineEdit(this);
  seedfilebtn = new QPushButton(QIcon(":/icons/browse.png"), " Browse", this);
  seedfilebtn->setFlat(true);
  seedfilebtn->setIconSize(QSize(20,20));

  // NACA 4 digits

  digitlbl = new QLabel("NACA digits", this);
  digitlbl->setEnabled(false);
  digitedit = new QLineEdit("0012");
  digitedit->setMaxLength(4);
  digitedit->setEnabled(false);

  // Top grid layout (basic settings)
  
  grid1 = new QGridLayout();
  grid1->addWidget(caselbl, 0, 0);
  grid1->addWidget(caseedit, 0, 1);
  grid1->addWidget(searchlbl, 1, 0);
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
  grid1->addWidget(seedlbl, 4, 0);
  grid1->addWidget(seedbox, 4, 1);
  hbox3 = new QHBoxLayout();
  hbox3->addWidget(seedfilelbl);
  hbox3->addWidget(seedfilebox);
  grid1->addLayout(hbox3, 4, 2);
  grid1->addWidget(seedfilebtn, 4, 3);
  hbox4 = new QHBoxLayout();
  hbox4->addWidget(digitlbl);
  hbox4->addWidget(digitedit);
  grid1->addLayout(hbox4, 5, 2);

  // Advanced settings label

  advlbl = new QLabel("<b>Advanced settings</b>", this);

  // Shape functions

  shapelbl = new QLabel("Shape functions", this);
  shapebox = new QComboBox(this);
  shapebox->addItem("Hicks-Henne");
  shapebox->addItem("NACA");

  // Number of shape functions

  nshapetlbl = new QLabel("Number of functions for top surface", this);
  nshapeblbl = new QLabel("Number of functions for bottom surface", this);
  nshapetbox = new QSpinBox(this);
  nshapebbox = new QSpinBox(this);
  nshapetbox->setMinimum(1);
  nshapetbox->setValue(4);
  nshapebbox->setMinimum(1);
  nshapebbox->setValue(4);

  // Initial perturbation

  initperturblbl = new QLabel("Initial perturbation", this);
  initperturbedit = new QLineEdit(this);
  initperturbedit->setText("0.025");
  initperturbedit->setToolTip("Max change in airfoil surface " +
                              QString("during initialization"));

  // Restart write frequency

  restlbl = new QLabel("Restart write frequency", this);
  restbox = new QSpinBox(this);
  restbox->setMinimum(1);
  restbox->setValue(20);

  // Whether to write design progression

  writedesignbox = new QCheckBox("Write design progression to files", this);
  writedesignbox->setChecked(true);

  // Minimum bump width

  minbumplbl = new QLabel("Minimum bump width", this);
  minbumpbox = new QLineEdit(this);
  minbumpbox->setText("0.1");

  // Bottom grid layout (advanced settings)

  grid2 = new QGridLayout();
  grid2->addWidget(shapelbl, 0, 0);
  grid2->addWidget(shapebox, 0, 1);
  grid2->addWidget(nshapetlbl, 1, 0);
  grid2->addWidget(nshapetbox, 1, 1);
  grid2->addWidget(nshapeblbl, 2, 0);
  grid2->addWidget(nshapebbox, 2, 1);
  grid2->addWidget(initperturblbl, 3, 0);
  grid2->addWidget(initperturbedit, 3, 1);
  grid2->addWidget(minbumplbl, 4, 0);
  grid2->addWidget(minbumpbox, 4, 1);
  grid2->addWidget(restlbl, 5, 0);
  grid2->addWidget(restbox, 5, 1);
  grid2->addWidget(writedesignbox, 6, 0);

  // Box layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addLayout(grid1);
  vbox->addWidget(advlbl);
  vbox->addLayout(grid2);
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
    seedfilelbl->setEnabled(true);
    seedfilebox->setEnabled(true);
    seedfilebtn->setEnabled(true);
  }
  else if (idx == 1)
  {
    seedfilelbl->setEnabled(false);
    seedfilebox->setEnabled(false);
    seedfilebtn->setEnabled(false);
    digitlbl->setEnabled(true);
    digitedit->setEnabled(true);
  }
}
