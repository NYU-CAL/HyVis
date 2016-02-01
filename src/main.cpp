#include "controlwindow.h"
#include "viewer.h"
#include "querywindow.h"
#include "plot2dviewer.h"
#include <QApplication>

// NOTE: VERSION NUMBER IS DEFINED WITHIN
// CONTROL WINDOW CODE :  controlwindow.h

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    Config c;
    c.hide();

    QueryWindow q;
    q.hide();

    Plot2DViewer plt;
    plt.hide();

    Viewer v;
    v.c = &c;
    v.qw = &q;
    v.plt = &plt;
    v.setFocusPolicy(Qt::ClickFocus);
    v.setFocus(Qt::OtherFocusReason);
    v.show();

    ControlWindow w;
    w.config = &c;
    w.viewer = &v;
    w.plt = &plt;
    w.qw = &q;
    w.show();

    // Set up for signal passing
    v.ctrlwin = &w;
    w.updateVariableBtns();

    //w.loadConfigFile("/Users/georgewong/Sites/CAL/default.config");

    return a.exec();
}
