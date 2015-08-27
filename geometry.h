#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>

#include <QVector2D>
#include <QVector3D>
#include <QGLFunctions>
#include <QGLShaderProgram>

#include "colormaps.h"

class Geometry : protected QGLFunctions
{
public:
    Geometry();
    ~Geometry();

    void init();
    void drawGeometryLeft(QGLShaderProgram *program);
    void drawGeometryRight(QGLShaderProgram *program);
    void drawGrid(QGLShaderProgram *program);

    void loadGeometry(double *t_jph, double **r_iph, int Nt, int *Nr, int Nc);
    void setValue(bool isLeft, int Nt, int *Nr, int Nc, double ***cells, int q);
    void setCmapMinmax(double **minmax, int numq);
    void updateCmapBound(int n, double val);
    void setCmapStats(double gamma, double center, double slope);
    void cycleCmap();

    void getCmapValue(double val, QVector3D *color);

private:
    // vbo[0]  :  points
    // vbo[1]  :  colors (left/top)
    // vbo[2]  :  colors (right/bottom)
    // vbo[3]  :  pt. order (faces)
    // vbo[4]  :  pt. order (grid)
    GLuint vbos[5];
    int Npts, Ngpts;
    double gamma, center, slope;
    double *cmap_minmax;
    int colormap;

};

#endif // GEOMETRY_H
