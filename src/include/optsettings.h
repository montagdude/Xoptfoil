#pragma once

#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>

/******************************************************************************/
//
// Header for optimization settings class
//
/******************************************************************************/
class OptSettings : public QWidget
{
  private:

    QLineEdit *caseedit;
    QComboBox *searchbox;
    QComboBox *globalbox;
    QPushButton *globalbtn;
    QComboBox *localbox;
    QPushButton *localbtn;

  public:

    // Constructor

    OptSettings ( QWidget *parent = 0 );
};
