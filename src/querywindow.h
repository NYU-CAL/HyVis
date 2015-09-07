#ifndef QUERYWINDOW_H
#define QUERYWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>

#include <QDebug>

#include "math.h"

namespace Ui {
class QueryWindow;
}

class QueryWindow : public QWidget
{
    Q_OBJECT

public:
    explicit QueryWindow(QWidget *parent = 0);
    ~QueryWindow();

    void setPosition(double p, double q, double r);
    void updateNumValues(int n, std::vector<QString> names);
    void setValues(double *values);

private:
    Ui::QueryWindow *ui;
    double p,q,r;
    int numValues;
    double *values;
    int numWdgtsDisplayed;
    std::vector <QWidget *>wdgts;
    std::vector <QVBoxLayout *>vboxs;
    std::vector <QLineEdit *>ledts;
    std::vector <QLabel *>lbls;
};

#endif // QUERYWINDOW_H
