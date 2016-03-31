#include "settingsbrowser.h"
#include "settingswindow.h"

#include <QWidget>
#include <QSplitter>
#include <QScrollArea>
#include <QList>

/******************************************************************************/
//
// Constructor
//
/******************************************************************************/
SettingsWindow::SettingsWindow(QWidget *parent) : QSplitter(parent)
{
  QList<int> default_sizes;

  // Widgets contained in splitter
  
  settingsbrowser = new SettingsBrowser(this);
  scrollarea = new QScrollArea(this);
  setOrientation(Qt::Horizontal);
  addWidget(settingsbrowser);
  addWidget(scrollarea);

  // Default sizes for widgets
  
  default_sizes.append(220);
  default_sizes.append(680);
  setSizes(default_sizes);
}
