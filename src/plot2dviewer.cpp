#include "plot2dviewer.h"
#include "ui_plot2dviewer.h"

Plot2DViewer::Plot2DViewer(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Plot2DViewer)
{
    ui->setupUi(this);
    this->ui->plotwindow->addGraph();
    this->ui->plotwindow->plotLayout()->insertRow(0);
}

Plot2DViewer::~Plot2DViewer()
{
    delete ui;
}

void Plot2DViewer::updatePlot(std::vector<double> *x, std::vector<double> *y, bool isFirst, double y0, double y1, std::string var)
{
    int i;
    QVector<double> xx(x->size());
    QVector<double> yy(x->size());
    double xmin, xmax;
    xmin = x->at(0); xmax = x->at(0);

    for (i=0; i<(int)x->size(); ++i) {
        xx[i] = x->at(i);
        yy[i] = y->at(i);
        if (xx[i] < xmin) xmin = xx[i];
        if (xx[i] > xmax) xmax = xx[i];
    }

    this->ui->plotwindow->graph(0)->setData(xx,yy);
    this->ui->plotwindow->replot();

    if (isFirst || !isFirst) {
        this->ui->plotwindow->xAxis->setRange(xmin, xmax);
        this->ui->plotwindow->yAxis->setRange(y0, y1);
    }

    this->setTitle(var);
}

void Plot2DViewer::setTitle(std::string title)
{
    try { this->ui->plotwindow->plotLayout()->removeAt(0); }
    catch(int e) { ; }
    this->ui->plotwindow->plotLayout()->addElement(0,0, new QCPPlotTitle(this->ui->plotwindow, title.c_str()));
}
