#pragma once

#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QLabel>

/******************************************************************************/
//
// Header for optimization settings class
//
/******************************************************************************/
class OptSettings : public QWidget
{
  Q_OBJECT

  private:

    QLineEdit *caseedit;
    QComboBox *searchbox;
    QComboBox *globalbox;
    QPushButton *globalbtn;
    QComboBox *localbox;
    QPushButton *localbtn;
    QComboBox *seedbox;
    QLabel *seedfilelbl;
    QLineEdit *seedfilebox;
    QPushButton *seedfilebtn;
    QLabel *digitlbl;
    QLineEdit *digitedit;
    QComboBox *shapebox;
    QSpinBox *nshapetbox;
    QSpinBox *nshapebbox;
    QLineEdit *initperturbedit;
    QLineEdit *minbumpbox;
    QSpinBox *restbox;
    QCheckBox *writedesignbox;

  private slots:

    void seedBoxChanged ( int idx );

  public:

    // Constructor

    OptSettings ( QWidget *parent = 0 );
};
