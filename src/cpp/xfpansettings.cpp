#include "xfpansettings.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QWidget>

/******************************************************************************/
//
// Constructor for Xfoil paneling settings class
//
/******************************************************************************/
XfPanSettings::XfPanSettings(QWidget *parent) : QWidget(parent)
{
  QLabel *lbl;
  QVBoxLayout *vbox;

  // Label for Xfoil paneling settings

  lbl = new QLabel("<b>Xfoil paneling settings</b>", this);

  // Layout

  vbox = new QVBoxLayout();
  vbox->addWidget(lbl);
  vbox->addStretch(0);

  setLayout(vbox);
}
