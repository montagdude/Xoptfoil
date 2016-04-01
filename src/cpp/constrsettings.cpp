#include "constrsettings.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for constraint settings class
//
/******************************************************************************/
ConstrSettings::ConstrSettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl;
  QVBoxLayout *vbox;

  // Label for constraints settings

  lbl = new QLabel("<b>Constraints settings</b>", this);

  // Layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addStretch(0);

  setLayout(vbox);
}
