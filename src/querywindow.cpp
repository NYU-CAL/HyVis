#include "querywindow.h"
#include "ui_querywindow.h"

QueryWindow::QueryWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::QueryWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("Query");
    this->numWdgtsDisplayed = 0;

    // Set up all possible value boxes
    int i;
    for (i=0; i<10; ++i) {
        QWidget *wdgt = new QWidget();
        QVBoxLayout *vbox = new QVBoxLayout();
        QLineEdit *ledt = new QLineEdit();
        QLabel *lbl = new QLabel("ace");

        ledt->setEnabled(false);
        vbox->addWidget(lbl);
        vbox->addWidget(ledt);
        wdgt->setLayout(vbox);

        this->wdgts.push_back(wdgt);
        this->vboxs.push_back(vbox);
        this->ledts.push_back(ledt);
        this->lbls.push_back(lbl);
        this->ui->horValLayout->addWidget(wdgt);
        this->numWdgtsDisplayed++;
    }


}

QueryWindow::~QueryWindow()
{
    // TODO
    delete ui;
}

void QueryWindow::setPosition(double p, double q, double r)
{
    this->p = p;
    this->q = q;
    this->r = r;
    QString P,Q,R;
    P.sprintf("%6.2f", p);
    Q.sprintf("%6.2f", q);
    R.sprintf("%6.2f", r);
    this->ui->editP->setText( P );
    this->ui->editQ->setText( Q );
    this->ui->editR->setText( R );
}

void QueryWindow::updateNumValues(int n, std::vector<QString> names)
{
    int i;
    // Remove old elements
    for (i=0; i<this->numWdgtsDisplayed; ++i) {
        this->ui->horValLayout->removeWidget(this->wdgts.at(i));
    }

    // Add "new" elements
    this->numValues = n;
    for (i=0; i<n; ++i) {
        this->ui->horValLayout->addWidget(this->wdgts.at(i));
        this->lbls.at(i)->setText(names.at(i));
    }
}

void QueryWindow::setValues(double *values)
{
    // Maybe should do data validation here, but I think
    // there's no harm if I don't...
    int i;
    for (i=0; i<this->numValues; ++i) {
        QString v;
        v.sprintf("%6.2f", values[i]);
        this->ledts.at(i)->setText(v);
    }
}
