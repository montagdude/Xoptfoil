#pragma once

#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QLabel>

// Forward declarations

class SettingsWindow;

/******************************************************************************/
//
// Header for optimization settings class
//
/******************************************************************************/
class OptSettings : public QWidget
{
  Q_OBJECT

  private:

    // Stored pointer to settings window

    SettingsWindow *settingswindow;

    // Widgets
    
    QLineEdit *caseedit;
    QComboBox *searchbox;
    QLabel *globallbl;
    QComboBox *globalbox;
    QPushButton *globalbtn;
    QLabel *locallbl;
    QComboBox *localbox;
    QPushButton *localbtn;
    QComboBox *seedbox;
    QLabel *seedfilelbl;
    QLineEdit *seedfilebox;
    QPushButton *seedfilebtn;
    QLabel *digitlbl;
    QSpinBox *digitbox1, *digitbox2, *digitbox3, *digitbox4;
    QComboBox *shapebox;
    QSpinBox *nshapetbox;
    QSpinBox *nshapebbox;
    QDoubleSpinBox *initperturbbox;
    QDoubleSpinBox *minbumpbox;
    QSpinBox *restbox;
    QCheckBox *writedesignbox;

  private slots:

    void searchBoxChanged ( int idx );
    void globalBoxChanged ( int idx );
    void globalBtnClicked ();
    void localBtnClicked ();
    void seedBoxChanged ( int idx );

  public:

    // Constructor

    OptSettings ( SettingsWindow *parent = 0 );
};
