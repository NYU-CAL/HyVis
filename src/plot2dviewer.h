#ifndef PLOT2DVIEWER_H
#define PLOT2DVIEWER_H

#include <QWidget>

namespace Ui {
class Plot2DViewer;
}

class Plot2DViewer : public QWidget
{
    Q_OBJECT

public:
    explicit Plot2DViewer(QWidget *parent = 0);
    ~Plot2DViewer();
    void updatePlot(std::vector<double> *x, std::vector<double> *y, bool isFirst, double y0, double y1, std::string var);
    void setTitle(std::string title);

private:
    Ui::Plot2DViewer *ui;
};

#endif // PLOT2DVIEWER_H
