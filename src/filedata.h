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


#ifndef FILEDATA_H
#define FILEDATA_H

#define FILEDATA_VERBOSE 1

#include <hdf5.h>
#include <cstdio>

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>

class FileData
{
public:
    FileData();
    ~FileData();

    int get_nc();
    int get_nq();
    double **get_minmax();
    double **get_cells();
    double *get_bounds();
    bool fileLoaded;
    
    virtual bool loadFromFile(const char *filename)=0;
    virtual void getValuesAt(double x, double y, double *values)=0;
    virtual int getValuesAtTheta(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)=0;
    virtual int getValuesAtR(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)=0;

    virtual void genGridData(double **gp, int *ngp, int **ci, int *nci, int **gpi, int *ngpi)=0;

protected:

    int numCells, numQ;
    double time;
    double **cells; // [i][j] -> i'th value of j'th cell in grid
    double **minmax; // [ [min0, min1, min2, ...] , [max0, max1, max2, ...] ] <-> [type,q]
    double bounds[4]; // [xmin xmax ymin ymax]
    
    void readPatch(hid_t *h5dst , void *data, hid_t type, int * start, int *loc_size, int *glo_size);

private:

    void freeAll();

};

#endif // FILEDATA_H
