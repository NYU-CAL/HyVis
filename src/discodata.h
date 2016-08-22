/*
 * We define the JET code grid (h5file format) below:
 *
 *  The root of the file itself has two parts:
 *      1. Data
 *      2. Grid
 *
 *  The Data section has one key: Cells. A two-
 *  dimensional array is stored at FILE[Data,Cells],
 *  which contains
 *
 *  The Grid section contains information about the
 *  general format of the instance of the grid. Its
 *  attributes are:
 *      1. Index
 *      2. Nr
 *      3. T
 *      4. p_kph
 *      5. t_jph
 *
 *  FILE[Grid,Index]
 *
 *  FILE[Grid,Nr]
 *
 *  FILE[Grid,T] corresponds to the time (in code units)
 *  of the h5 file.
 *
 *  FILE[Grid,p_kph] gives the phi values bounding the cells
 *
 *  FILE[Grid,t_jph] gives the theta values corresponding
 *  to the "edges" of the cells
 *
 *
 */


#ifndef DISCODATA_H
#define DISCODATA_H

#include <hdf5.h>
#include <cstdio>

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>

class DiscoData : public FileData
{
public:
    DiscoData();
    ~DiscoData();

    bool loadFromFile(const char *filename);
    void getValuesAt(double x, double y, double *values);
    int getValuesAtTheta(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata);
    int getValuesAtR(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata);

    void genGridData(QVector2D **gp, int *ngp, GLuint **ci, int *nci, 
                        GLuint **gpi, int *ngpi);

private:

    // These all relate to DISCO specific variables
    int numR, numZ;
    int *Np;         // An array of numR #-of-phi-zones values
    double *r_jph;   // An array of numR+1 radial values
    double *z_kph;   // An array of numZ+1 vertical values
    double **p_iph;  // [j][i] -> i'th azimuthal position of j'th annulus
    void freeAllDisco();
};

#endif // FILEDATA_H
