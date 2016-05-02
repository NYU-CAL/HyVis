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

#define FILEDATA_VERBOSE 0

#include <hdf5.h>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>
//#include "unit.h"

class Unit{
 public:
  Unit();
  ~Unit();
  int nums;
  double M;
  double L;
  double T;
  double Rho;
  double Pre;
  double V;
  double *vars;
};


class FileData
{
public:
    FileData();
    ~FileData();

    bool loadFromFile(const char *filename);
    bool ToCGS();    
    int get_nc();
    int get_nt();
    int get_nq();
    int *get_nr();
    double *get_tjph();
    double **get_riph();
    double **get_minmax();
    double ***get_cells();
    void getValuesAt(double x, double y, double *values);
    int getValuesAtTheta(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata);
    int getValuesAtR(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata);

    bool fileLoaded;
    
private:
    
    // These all relate to JET specific variables
    int numPhi, numTheta;
    int numCells, numQ;
    double time;
    int *Nr;         // An array of numTheta #-of-radial-zones values
    double *t_jph;   // An array of numTheta+1 angle values
    double **r_iph;  // [i][j] -> j'th radial position of i'th ray
    double ***cells; // [i][j][k] -> i'th value of k'th cell in j'th ray
    double **minmax; // [ [min0, min1, min2, ...] , [max0, max1, max2, ...] ] <-> [type,q]

    void readPatch(hid_t *h5dst , void *data, hid_t type, int * start, int *loc_size, int *glo_size);

    void freeAll();
    Unit unit;
};

#endif // FILEDATA_H
