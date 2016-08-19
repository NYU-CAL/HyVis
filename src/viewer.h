#ifndef VIEWER_H
#define VIEWER_H

#include "geometry.h"
#include "filedata.h"
#include "jetdata.h"
#include "querywindow.h"
#include "plot2dviewer.h"
#include <controlwindow.h>
#include <fstream>

#include <cstdio>
#include <cstdlib>

#include <QGLWidget>
#include <QGLShaderProgram>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>

class Config;
class ControlWindow;

const double DEG_RAD = M_PI / 180.;
const int INPUT_LENGTH = 17;
enum INPUT {MOUSE_LEFT, MOUSE_RIGHT, KEY_SHIFT, KEY_CTRL, KEY_LEFT, KEY_RIGHT,  // 6
            KEY_W, KEY_A, KEY_S, KEY_D, KEY_E, KEY_C, KEY_Q, KEY_G, KEY_MINUS,  // 9
            KEY_PLUS, KEY_L};            // 2

class Viewer : public QGLWidget
{
    Q_OBJECT
public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();

    Config *c;
    QueryWindow *qw;
    Plot2DViewer *plt;
    ControlWindow *ctrlwin;

    bool loadFile(const char *filename);
    void setByAxes(double x0, double x1, double y0, double y1);
    void setCenter(double x, double y);
    void setXScale(double xscale);
    void setYScale(double yscale);
    void squareScale();

    void setGridVisible(bool isLeft, bool value);
    void rotate(double angle);
    void setIsMirrored(bool mirror);
    void setIsDoubleMirrored(bool doublemirrored);
    void updateOrthoMatrix();
    void updateOrthoMatrix(bool updateControl);
    void setBgColor(double r, double g, double b);
    void setLeftValue(int n);
    void setRightValue(int n);
    int get_nq();
    void rescale(int i);

    // Colormap
    void showScale(bool value);
    void setLogscale(bool log);
    void showColorbar(bool value);
    void setGammaValue(double gamma);
    void setCenterValue(double center);
    void setSlopeValue(double slope);
    void setCmapMinmax(double min, double max);
    void requestCboundReset();
    void cycleCmap();

    void set1DPlot(bool visible, bool isTheta);
    void replot1DPlot();

    void savePPM(const char *fname);

    int leftDisplayVar, rightDisplayVar;

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintEvent(QPaintEvent *event);

    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent* event);
    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

private:
    bool glDoneInit;
    Geometry g; // What did the acorn say when it grew up?
    JetData filedata;
    QGLShaderProgram program;
    QGLShaderProgram whiteprogram;
    QMatrix4x4 matrix;
    QMatrix4x4 matrixRight;
    QMatrix4x4 invmatrix;
    QMatrix4x4 invmatrixRight;

    bool input[INPUT_LENGTH];
    double pointCenter[2];
    double pointScale[2];
    double clickPoint[2];
    double clickCenter[2];
    int bgColor[3];

    bool drawBox;
    bool drawScale;
    bool drawColorbar;
    bool drawGridL, drawGridR;
    double drawBoxCoordinates[4];
    double *valuesAtPointer;

    // Visualization parameters
    bool isMirrored, isDoubleMirrored;
    double rotateAngle;
    float gridColor[3];

    // 1-D Plotting
    void update1DPlot(int x, int y, bool isFirst);
    bool drawPlotLine, plotIsTheta;
    double plotLineData[3];

    void initShaders();

    // Note that this function does not account for pointCenter ... you
    // have to add that manually! This is so that the click / drag feature work.
    void pxToPoint(int x, int y, double *tx, double *ty);
    void screenToData(int sx, int sy, double *tx, double *ty);


signals:

public slots:

};

#endif // VIEWER_H
