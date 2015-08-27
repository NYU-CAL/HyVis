#include "config.h"
#include "ui_config.h"

Config::Config(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Config)
{
    ui->setupUi(this);
}

Config::~Config()
{
    delete ui;
}

void Config::setNq(int nq)
{
    if (nq > 0) this->ui->var0form->setDisabled(true);
    if (nq > 1) this->ui->var1form->setDisabled(true);
    if (nq > 2) this->ui->var2form->setDisabled(true);
    if (nq > 3) this->ui->var3form->setDisabled(true);
    if (nq > 4) this->ui->var4form->setDisabled(true);
    if (nq > 5) this->ui->var5form->setDisabled(true);
    if (nq > 6) this->ui->var6form->setDisabled(true);
}

QString Config::getName(int n)
{
    if (n == 0) return this->ui->var0name->text();
    if (n == 1) return this->ui->var1name->text();
    if (n == 2) return this->ui->var2name->text();
    if (n == 3) return this->ui->var3name->text();
    if (n == 4) return this->ui->var4name->text();
    return QString("");
}
