#include "plot2dviewer.h"
#include "ui_plot2dviewer.h"

Plot2DViewer::Plot2DViewer(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Plot2DViewer)
{
    ui->setupUi(this);
    this->ui->plotwindow->addGraph();
}

Plot2DViewer::~Plot2DViewer()
{
    delete ui;
}

void Plot2DViewer::updatePlot(std::vector<double> *x, std::vector<double> *y, bool isFirst)
{
    int i;
    QVector<double> xx(x->size());
    QVector<double> yy(x->size());
    double xmin, xmax, ymin, ymax;
    xmin = x->at(0); xmax = x->at(0);
    ymin = y->at(0); ymax = y->at(0);

    for (i=0; i<(int)x->size(); ++i) {
        xx[i] = x->at(i);
        yy[i] = y->at(i);
        if (xx[i] < xmin) xmin = xx[i];
        if (xx[i] > xmax) xmax = xx[i];
        if (yy[i] < ymin) ymin = yy[i];
        if (yy[i] > ymax) ymax = yy[i];
    }

    this->ui->plotwindow->graph(0)->setData(xx,yy);
    this->ui->plotwindow->replot();

    if (isFirst || !isFirst) {
        this->ui->plotwindow->xAxis->setRange(xmin, xmax);
        this->ui->plotwindow->yAxis->setRange(0.0, 16.0); //ymin, ymax);
    }


    /*
    // give the axes some labels:
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");
    */
}
