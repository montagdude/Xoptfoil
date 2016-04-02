#pragma once

#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
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
    QPushButton *seedfilebtn;
    QLabel *digitlbl;
    QLineEdit *digitedit;

  private slots:

    void seedBoxChanged ( int idx );

  public:

    // Constructor

    OptSettings ( QWidget *parent = 0 );
};
