#include "xfanasettings.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for Xfoil analysis settings class
//
/******************************************************************************/
XfAnaSettings::XfAnaSettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl;
  QVBoxLayout *vbox;

  // Label for Xfoil analysis settings

  lbl = new QLabel("<b>Xfoil analysis settings</b>", this);

  // Layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addStretch(0);

  setLayout(vbox);
}
