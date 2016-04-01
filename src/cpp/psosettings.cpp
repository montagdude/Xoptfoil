#include "psosettings.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for particle swarm settings class
//
/******************************************************************************/
PSOSettings::PSOSettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl;
  QVBoxLayout *vbox;

  // Label for particle swarm settings

  lbl = new QLabel("<b>Particle swarm settings</b>", this);

  // Layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addStretch(0);

  setLayout(vbox);
}
