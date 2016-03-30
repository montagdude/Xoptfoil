#pragma once

#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>

/******************************************************************************/
//
// Header for settings browser class
//
/******************************************************************************/
class SettingsBrowser : public QListWidget
{
  private:

    // List items
    
    QListWidgetItem *optimitem;
    QListWidgetItem *operitem;
    QListWidgetItem *constritem;
    QListWidgetItem *inititem;
    QListWidgetItem *psoitem;
    QListWidgetItem *gaitem;
    QListWidgetItem *simplexitem;
    QListWidgetItem *xfrunitem;
    QListWidgetItem *xfpanitem;

  public:

    // Constructor

    SettingsBrowser ( QWidget *parent = 0 );
};
