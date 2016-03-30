#pragma once

#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>

/******************************************************************************/
//
// Header for navigator list class
//
/******************************************************************************/
class NavigatorList : public QListWidget
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

    NavigatorList ( QWidget *parent = 0 );
};
