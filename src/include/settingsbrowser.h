#pragma once

#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>

// Forward declarations

class SettingsWindow;

/******************************************************************************/
//
// Header for settings browser class
//
/******************************************************************************/
class SettingsBrowser : public QListWidget
{
  private:

    // Stored pointer to settings window

    SettingsWindow *settingswindow;

    // List items
    
    QListWidgetItem *optimitem;
    QListWidgetItem *operitem;
    QListWidgetItem *constritem;
    QListWidgetItem *inititem;
    QListWidgetItem *psoitem;
    QListWidgetItem *gaitem;
    QListWidgetItem *simplexitem;
    QListWidgetItem *xfanaitem;
    QListWidgetItem *xfpanitem;

  private slots:

    // Occurs when selection is changed

    void selectionChanged ();

  public:

    // Constructor

    SettingsBrowser ( SettingsWindow *parent = 0 );
};
