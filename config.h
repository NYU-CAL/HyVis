#ifndef CONFIG_H
#define CONFIG_H

#include <QDialog>

namespace Ui {
class Config;
}

class Config : public QDialog
{
    Q_OBJECT

public:
    explicit Config(QWidget *parent = 0);
    ~Config();
    void setNq(int nq);
    QString getName(int n);

private:
    Ui::Config *ui;
};

#endif // CONFIG_H
