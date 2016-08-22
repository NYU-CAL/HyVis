#ifndef GEOMETRY_H
#define GEOMETRY_H

#define GEOMETRY_VERBOSE 0

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

    void setLogscale(bool log);
    void loadGeometry(QVector2D *gp, int ngp, int nc, GLuint *ci, int nci, 
                        GLuint *gpi, int ngpi);
    void setValue(bool isLeft, int Nc, double **cells, int q);
    void setCmapMinmax(double **minmax, int numq);
    void updateCmapBound(int n, double val);
    void setCmapStats(double gamma, double center, double slope);
    void cycleCmap();

    void getCmapValue(double val, QVector3D *color);
    void unitToTrue(double in, double *out);

private:
    // vbo[0]  :  points
    // vbo[1]  :  colors (left/top)
    // vbo[2]  :  colors (right/bottom)
    // vbo[3]  :  pt. order (faces)
    // vbo[4]  :  pt. order (grid)
    GLuint vbos[5];
    int Ndpts, Ngpts, Nglns, Ndpts_tot;
    double gamma, center, slope;
    double *cmap_minmax;
    int colormap;
    int leftQ;
    bool log;

};

#endif // GEOMETRY_H
