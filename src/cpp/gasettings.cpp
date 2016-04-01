#include "gasettings.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for genetic algorithm settings class
//
/******************************************************************************/
GASettings::GASettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl;
  QVBoxLayout *vbox;

  // Label for genetic algorithn settings

  lbl = new QLabel("<b>Genetic algorithm settings</b>", this);

  // Layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addStretch(0);

  setLayout(vbox);
}
