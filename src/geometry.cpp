#include "geometry.h"

Geometry::Geometry()
{
    //this->Npts = 0;
    this->Ngpts = 0;
    this->Nglns = 0;
    this->Ndpts = 0;
    this->Ndpts_tot = 0;
    this->gamma = 0.8;
    this->center = 0.5;
    this->slope = 1.0;
    this->cmap_minmax = NULL;
    this->colormap = 2;
    this->leftQ = -1;
}

Geometry::~Geometry()
{
    glDeleteBuffers(4, this->vbos);
    // TODO
}

void Geometry::init()
{
    initializeGLFunctions(); //TODO: Is this necessary?
    glGenBuffers(5, this->vbos);

    // For now, because we assume there's no geometry,
    // we initialize a square.
    QVector2D vertices[] = {
        QVector2D( 0.0f, -1.0f),
        QVector2D(-1.0f, -1.0f),
        QVector2D(-1.0f,  1.0f),
        QVector2D( 0.0f,  1.0f),

        QVector2D( 1.0f, -1.0f),
        QVector2D( 0.0f, -1.0f),
        QVector2D( 0.0f,  1.0f),
        QVector2D( 1.0f,  1.0f),
    };
    QVector3D colors[] = {
        QVector3D( 1.0f, 0.0f, 0.0f ),
        QVector3D( 1.0f, 0.0f, 0.0f ),
        QVector3D( 1.0f, 0.0f, 0.0f ),
        QVector3D( 1.0f, 0.0f, 0.0f ),

        QVector3D( 0.0f, 0.0f, 1.0f ),
        QVector3D( 0.0f, 0.0f, 1.0f ),
        QVector3D( 0.0f, 0.0f, 1.0f ),
        QVector3D( 0.0f, 0.0f, 1.0f ),
    };

    GLuint indices[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

    this->Ngpts = 8;
    this->Nglns = 8;
    this->Ndpts = 8;
    this->Ndpts_tot = 8;

    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[0]);
    glBufferData(GL_ARRAY_BUFFER, this->Ngpts * sizeof(QVector2D), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[1]);
    glBufferData(GL_ARRAY_BUFFER, this->Ndpts_tot * sizeof(QVector3D), colors, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[2]);
    glBufferData(GL_ARRAY_BUFFER, this->Ndpts_tot * sizeof(QVector3D), colors, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbos[3]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->Ndpts * sizeof(GLuint), indices, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Geometry::drawGeometryLeft(QGLShaderProgram *program)
{
    int vertexLocation = program->attributeLocation("a_position");
    int colorLocation = program->attributeLocation("a_color");

    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[0]);
    program->enableAttributeArray(vertexLocation);
    glVertexAttribPointer(vertexLocation, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), NULL);

    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[1]);
    program->enableAttributeArray(colorLocation);
    glVertexAttribPointer(colorLocation, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), NULL);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbos[3]);

    glDrawElements(GL_QUADS, this->Ndpts, GL_UNSIGNED_INT, 0);
}

void Geometry::drawGeometryRight(QGLShaderProgram *program)
{
    int vertexLocation = program->attributeLocation("a_position");
    int colorLocation = program->attributeLocation("a_color");

    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[0]);
    program->enableAttributeArray(vertexLocation);
    glVertexAttribPointer(vertexLocation, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), NULL);

    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[2]);
    program->enableAttributeArray(colorLocation);
    glVertexAttribPointer(colorLocation, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), NULL);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbos[3]);

    glDrawElements(GL_QUADS, this->Ndpts, GL_UNSIGNED_INT, 0);

}

void Geometry::drawGrid(QGLShaderProgram *program)
{
    int vertexLocation = program->attributeLocation("a_position");
    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[0]);
    program->enableAttributeArray(vertexLocation);
    glVertexAttribPointer(vertexLocation, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), NULL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbos[4]);
    glDrawElements(GL_LINES, this->Nglns, GL_UNSIGNED_INT, 0);
}

void Geometry::setLogscale(bool log)
{
    this->log = log;
}

//void Geometry::loadGeometry(double *t_jph, double **r_iph, int Nt, int *Nr, int Nc)
void Geometry::loadGeometry(double *gp, int ngp, int nc, int *ci, int nci, 
                            int *gpi, int ngpi)
{
// gp: array of (x,y) locations of grid vertices: xi=gp[2*i], yi=[2*i+1]
// ngp: number of grid vertices. length of gp is 2*ngp
// nc: Total number of grid cells
// ci: Array of indices (for "cells") to draw. 4 PER CELL.
// nci: Number of cell indices (length of ci)
// gpi: Grid Point Indices, indices of gp to draw
// ngpi: length of gpi, equals twice the number of grid lines drawn.

    int i,j, c,d;

    this->Ngpts = ngp;
    this->Nglns = ngpi;
    this->Ndpts_tot = 4*nc;
    this->Ndpts = nci;

    QVector2D *vertices = (QVector2D *)malloc(sizeof(QVector2D) * this->Ngpts );
    QVector3D   *colors = (QVector3D *)malloc(sizeof(QVector3D) * this->Ndpts_tot );
    GLuint     *indices = (GLuint *)   malloc(sizeof(GLuint)    * this->Ndpts );
    GLuint    *gindices = (GLuint *)   malloc(sizeof(GLuint)    * this->Nglns);
    c = 0;
    for(c=0; c<this->Ngpts; c++)
        vertices[c] = QVector2D((float)gp[2*c], (float)gp[2*c+1]);
    for(c=0; c<nci; c++)
        indices[c] = ci[c];
    for(c=0; c<this->Ndpts_tot; c++)
        colors[c] = QVector3D(1.0f, 1.0f, 1.0f);
    for(c=0; c<this->Nglns; c++)
        gindices[c] = gpi[c];

    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[0]);
    glBufferData(GL_ARRAY_BUFFER, this->Ngpts * sizeof(QVector2D), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[1]);
    glBufferData(GL_ARRAY_BUFFER, this->Ndpts_tot * sizeof(QVector3D), colors, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbos[2]);
    glBufferData(GL_ARRAY_BUFFER, this->Ndpts_tot * sizeof(QVector3D), colors, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbos[3]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->Ndpts * sizeof(GLuint), indices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbos[4]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->Nglns * sizeof(GLuint), gindices, GL_STATIC_DRAW);

    // Housekeeping
    free(vertices);
    free(colors);
    free(indices);
    free(gindices);
}

void Geometry::setCmapStats(double gamma, double center, double slope)
{
    if (gamma < 0.0)  ; else this->gamma = gamma;
    if (center < 0.0) ; else this->center = center;
    if (slope < 0.0)  ; else this->slope = slope;
}

void Geometry::cycleCmap()
{
    this->colormap = (this->colormap + 1) % CMAP::NUM;
}

void Geometry::setValue(bool isLeft, int Nc, double **cells, int q)
{
    QVector3D *colors   = (QVector3D *)malloc(sizeof(QVector3D) * Nc * 4);

    double x0 = this->cmap_minmax[0+2*q];
    double x1 = this->cmap_minmax[1+2*q];

    if (this->log) {
        if (x0 < 1e-4) x0 = 1e-4; // This is a bit of a magic number...
        x0 = log10(x0);
        x1 = log10(x1);
    }

    double m = 0.0;
    if ( (x1-x0) < 0.0000000001 ) {
        x0 = x0 - 0.5;
        m = 1.0;
    }
    else m = 1.0 / (x1 - x0);

    //int i, j, c=0;
    int i, c=0;
    for (i=0; i<Nc; ++i) {
    //for (i=0; i<Nt; ++i) {
        //for (j=1; j<Nr[i]; ++j) {
            //double value = cells[q][i][j];
            double value = cells[q][i];

            if (this->log) value = log10(value);

            value -= x0;
            value *= m;

            if (value < 0) value = 0;
            if (value > 1) value = 1;

            QVector3D col;
            this->getCmapValue(value, &col);
            colors[c] = QVector3D(col.x(), col.y(), col.z()); c++;
            colors[c] = QVector3D(col.x(), col.y(), col.z()); c++;
            colors[c] = QVector3D(col.x(), col.y(), col.z()); c++;
            colors[c] = QVector3D(col.x(), col.y(), col.z()); c++;
        //}
    }

    if (isLeft) {
        glBindBuffer(GL_ARRAY_BUFFER, this->vbos[1]);
        this->leftQ = q;
    }
    else glBindBuffer(GL_ARRAY_BUFFER, this->vbos[2]);

    glBufferData(GL_ARRAY_BUFFER, this->Ndpts_tot * sizeof(QVector3D), colors, GL_STATIC_DRAW);

    free(colors);
}

void Geometry::setCmapMinmax(double **minmax, int numq)
{
    int i;
    if (this->cmap_minmax != NULL) free(this->cmap_minmax);
    this->cmap_minmax = (double *)malloc(numq * 2 * sizeof(double));
    for (i=0; i<numq; ++i) {
        this->cmap_minmax[0 + 2*i] = minmax[0][i];
        this->cmap_minmax[1 + 2*i] = minmax[1][i];
    }
}

void Geometry::updateCmapBound(int n, double val)
{
    if (this->cmap_minmax != NULL) this->cmap_minmax[n] = val;
}

void Geometry::getCmapValue(double val, QVector3D *color)
{
    // ////////////////////////////// //
    // Rescale according to transform //
    // ////////////////////////////// //
    // Compute a two-part spline using the hermite polynomials.
    // The matrix is:
    //  2 -2  1  1  f(0)
    // -3  3 -2 -1  f(1)
    //  0  0  1  0  f'(0)
    //  1  0  0  0  f'(1)

    double mid = 0.5;
    double a0,b0,c0,d0, a1,b1,c1,d1;

    a0 = -2.0 * mid + 1.0 / this->slope + this->slope;
    b0 = 3.0 * mid - 2.0 / this->slope - this->slope;
    c0 = 1.0 / this->slope; d0 = 0.0;

    a1 = 2.0 * mid - 2.0  + this->slope + 1.0 / this->slope;
    b1 = -3.0 * mid + 3.0 - 2.0 * this->slope - 1.0 / this->slope;
    c1 = this->slope; d1 = mid;

    if (val < this->center) {
        val /= this->center;
        val = a0*val*val*val + b0*val*val + c0*val + d0;
    } else {
        val -= this->center;
        val /= (1.0 - this->center);
        val = a1*val*val*val + b1*val*val + c1*val + d1;
    }

    // //////////////////// //
    // Perform colormapping //
    // //////////////////// //
    double r,g,b;

    // For reference:
    //  0 -> pseudo-jet
    //  1 -> hot
    //  2 -> OPTION_A
    //  3 -> OPTION_C
    //  4 -> OPTION_D

    switch (this->colormap) {

    case 0:
    {
        double nexp = 8.0;
        r = exp(-nexp*pow(val-5./6.,2.0)) + .25*exp(-nexp*pow(val+1./6.,2.0));
        g = exp(-nexp*pow(val-3./6.,2.0));
        b = exp(-nexp*pow(val-1./6.,2.0)) + .25*exp(-nexp*pow(val-7./6.,2.0));
        break;
    }

    case 1:
    {
        if (val < 0.4) { r = val / 0.4; g = 0.0; b = 0.0; }
        else if (val < 0.8) { g = (val - 0.4 ) / 0.4; r = 1.0; b = 0.0; }
        else { r = 1.0; g = 1.0; b = (val - 0.8) / 0.2; }
        break;
    }

    case 2:
    {
        // TODO, do interpolation between cmap values
        int idx = (int) ( 255 * val );
        if (idx < 0) idx = 0;
        if (idx > 255) idx = 255;
        r = CMAP::OPTION_A[idx][0];
        g = CMAP::OPTION_A[idx][1];
        b = CMAP::OPTION_A[idx][2];
        break;
    }

    case 3:
    {
        int idx = (int) ( 255 * val );
        r = CMAP::OPTION_C[idx][0];
        g = CMAP::OPTION_C[idx][1];
        b = CMAP::OPTION_C[idx][2];
        break;
    }

    case 4:
    {
        int idx = (int) ( 255 * val );
        r = CMAP::OPTION_D[idx][0];
        g = CMAP::OPTION_D[idx][1];
        b = CMAP::OPTION_D[idx][2];
        break;
    }

    default:
    {
        r = 0; g = 0; b = 0;
        break;
    }

    }

    // Perform gamma correction and set the color
    r = pow(r, this->gamma);
    g = pow(g, this->gamma);
    b = pow(b, this->gamma);

    if (r < 0) r = 0.0; if (r > 1) r = 1.0;
    if (g < 0) g = 0.0; if (g > 1) g = 1.0;
    if (b < 0) b = 0.0; if (b > 1) b = 1.0;

    color->setX( (float)r );
    color->setY( (float)g );
    color->setZ( (float)b );
}

void Geometry::unitToTrue(double in, double *out)
{
    double min, max;

    if (this->leftQ == -1) {
        min = 0.0;
        max = 1.0;
    } else {
        min = this->cmap_minmax[0+2*this->leftQ];
        max = this->cmap_minmax[1+2*this->leftQ];
    }

    if (this->log) {
        if (min < 1e-4) min = 1e-4;
        min = log10(min);
        max = log10(max);
    }

    *out = min + (max-min) * in;
    if (this->log) *out = pow(10.0, *out);
}
