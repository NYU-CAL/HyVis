#include "filedata.h"

Unit::Unit()
{
  M = 2e33;
  L = 6.8e10;
  T = 96.1665;
  Rho = 6.36067575819255;
  Pre = 3.1803378790962755e18;
  V = 7071.0678118654751;
  nums = 10; // note this number should larger than numQ;
  vars = (double *) malloc(sizeof(double)*nums);
  for(int i=0; i<nums; i++){
    vars[i] = 1.0;
  }
  vars[0] = Rho;
  vars[1] = Pre;
  vars[2] = V;
  vars[3] = V;
  vars[4] = V;
}

Unit::~Unit()
{
  free(vars);
}
    

FileData::FileData()
{
    this->numCells = 0;
    this->numPhi = 1;
    this->numTheta = 0;
    this->time = 0.0;
    this->fileLoaded = false;

    this->minmax = (double **)malloc(sizeof(double *) * 2);

    // Initialize null pointers;
    this->cells = NULL;
    this->r_iph = NULL;
    this->t_jph = NULL;
    this->Nr = NULL;
    this->minmax[0] = NULL;
    this->minmax[1] = NULL;
}

FileData::~FileData()
{
    this->freeAll();
    free(this->minmax);
}

bool FileData::ToCGS(){
  
  this->time = this->time*this->unit.T;
  if( this->numQ > this->unit.nums ) printf("Warning Unit Convert Variable number less than the actual variable number\n");
  //  printf("%.4e   %.4e   %.4e   %.4e\n", this->unit.L, this->unit.Rho, this->unit.Pre, this->unit.V);
  for (int i=0; i<this->numTheta; ++i) {
    for (int j=0; j<Nr[i]; ++j) {
      //      this->r_iph[i][j] = this->r_iph[i][j]*this->unit.L;
      for (int k=0; k<this->numQ; ++k) {
	this->cells[k][i][j] = this->cells[k][i][j]*this->unit.vars[k];
      }
    }
  }
  for(int k=0; k<this->numQ; ++k){
    this->minmax[0][k] = this->minmax[0][k]*this->unit.vars[k];
    this->minmax[1][k] = this->minmax[1][k]*this->unit.vars[k];
  }
  return true;
}


bool FileData::loadFromFile(const char *filename)
{
    // Useful variables
    int i, j;
    hsize_t dims[3];
    int start[2];
    int loc_size[2];
    int glo_size[2];
    int *tindex;

    #if FILEDATA_VERBOSE
    fprintf(stderr, "Loading data from file: %s\n", filename);
    #endif

    hid_t h5file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5file < 0) {
        fprintf(stderr, "Unable to open h5 file.\n");
        return false;
    }

    hid_t h5gridgroup = H5Gopen1(h5file, "Grid");
    hid_t h5datagroup = H5Gopen1(h5file, "Data");

    if (h5gridgroup < 0 || h5datagroup < 0) {
        fprintf(stderr, "Successfully opened file, but encountered unexpected format.\n");
        return false;
    }

    // For alpha purposes, we can assume this is the commit point -- and we'll free all memory...
    this->fileLoaded = false;
    this->freeAll();

    // First step is loading the grid information

        // Get time
        hid_t h5_t = H5Dopen1(h5gridgroup, "T");
        H5Dread(h5_t, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)&(this->time));
        H5Dclose(h5_t);

        // Get number of phi plates
        this->numPhi = 1;   // Technically, I could load this; however, becaues the code
                            // doesn't support nphi>1, I'll leave it as is.

        // Get number of theta rays
        hid_t h5_theta = H5Dopen1(h5gridgroup, "t_jph");
        hid_t h5_theta_spc = H5Dget_space(h5_theta);
        H5Sget_simple_extent_dims(h5_theta_spc, dims, NULL);
        this->numTheta = dims[0] - 1; // TODO ... setting things here causes issues.
        H5Sclose(h5_theta_spc);

        // Then we actually load the t_jph angles
        this->t_jph = (double *)malloc(dims[0] * sizeof(double));
        start[0] = 0; start[1] = 0;
        loc_size[0] = dims[0]; loc_size[1] = 1;
        glo_size[0] = dims[0]; glo_size[1] = 1;
        readPatch(&h5_theta, this->t_jph, H5T_NATIVE_DOUBLE, start, loc_size, glo_size);
        H5Dclose(h5_theta);

        // ... and we get the number of radial cells per theta ray
        this->Nr = (int *)malloc(this->numTheta * sizeof(int));
        hid_t h5_Nr = H5Dopen1(h5gridgroup, "Nr");
        loc_size[0]--; glo_size[0]--;
        readPatch(&h5_Nr, this->Nr, H5T_NATIVE_INT, start, loc_size, glo_size);
        H5Dclose(h5_Nr);

        // ... and finally, we load the indices for the cells, into a temp array.
        tindex = (int *)malloc(this->numTheta * sizeof(int));
        hid_t h5_tindex = H5Dopen1(h5gridgroup, "Index");
        readPatch(&h5_tindex, tindex, H5T_NATIVE_INT, start, loc_size, glo_size);
        H5Dclose(h5_tindex);

    // Then we load the information in the cells (geometry + values)
    hid_t h5_cells = H5Dopen1(h5datagroup, "Cells");
    hid_t h5_cells_spc = H5Dget_space(h5_cells);
    H5Sget_simple_extent_dims(h5_cells_spc, dims, NULL);
    this->numCells = dims[0];
    this->numQ = dims[1] - 1;
    H5Sclose(h5_cells_spc);

    // Set up minmax
    this->minmax[0] = (double *)malloc(sizeof(double) * this->numQ);
    this->minmax[1] = (double *)malloc(sizeof(double) * this->numQ);

    this->r_iph = (double **)malloc(sizeof(double *) * this->numTheta);
    this->cells = (double ***)malloc(sizeof(double **) * this->numQ);
    for (i=0; i<this->numQ; ++i) {
        cells[i] = (double **)malloc(sizeof(double *) * this->numTheta);
    }

    loc_size[1] = this->numQ+1;
    glo_size[0] = this->numCells; glo_size[1] = this->numQ + 1;
    start[0] = - this->Nr[0];
    for (i=0; i<this->numTheta; ++i) {
        this->r_iph[i] = (double *)malloc(sizeof(double) * Nr[i]);
        for (j=0; j<this->numQ; ++j) cells[j][i] = (double *)malloc(sizeof(double) * Nr[i]);
        loc_size[0] = this->Nr[i];
        start[0] = tindex[i];
        double buffer[this->Nr[i] * (this->numQ+1)];
        readPatch(&h5_cells, buffer, H5T_NATIVE_DOUBLE, start, loc_size, glo_size);

        // (this is annoying, but it's the best way to do it)
        if (i==0) {
            for (j=0; j<this->numQ; ++j) {
                this->minmax[0][j] = buffer[j];
                this->minmax[1][j] = buffer[j];
            }
        }

        // Load from the buffer into the proper class variables
        for (j=0; j<Nr[i]; ++j) {
            this->r_iph[i][j] = buffer[j * (this->numQ + 1) + this->numQ];
            for (int k=0; k<this->numQ; ++k) {
                this->cells[k][i][j] = buffer[j * (this->numQ + 1) + k];
                if (this->minmax[0][k] > this->cells[k][i][j]) this->minmax[0][k] = this->cells[k][i][j];
                if (this->minmax[1][k] < this->cells[k][i][j]) this->minmax[1][k] = this->cells[k][i][j];
            }
        }
    }

    H5Dclose(h5_cells);

    // Housekeeping!
    H5Gclose(h5gridgroup);
    H5Gclose(h5datagroup);
    H5Fclose(h5file);

    free(tindex);

    this->fileLoaded = true;

    #if FILEDATA_VERBOSE
    fprintf(stderr, "Done loading data.\n");
    #endif

    return true;
}

double *FileData::get_tjph()
{
    if (this->fileLoaded) return this->t_jph;
    return NULL;
}

double **FileData::get_riph()
{
    if (this->fileLoaded) return this->r_iph;
    return NULL;
}

double **FileData::get_minmax()
{
    if (this->fileLoaded) return this->minmax;
    return NULL;
}

int FileData::get_nc()
{
    if (this->fileLoaded) return this->numCells;
    return 0;
}

int FileData::get_nt()
{
    if (this->fileLoaded) return this->numTheta;
    return 0;
}

int FileData::get_nq()
{
    if (this->fileLoaded) return this->numQ;
    return 0;
}

int *FileData::get_nr()
{
    if (this->fileLoaded) return this->Nr;
    return NULL;
}

double ***FileData::get_cells()
{
    if (this->fileLoaded) return this->cells;
    return NULL;
}

// This function sets values[i] to the i'th q of the cell
// located at (x,y). If such a cell does not exist, it
// simply returns.
void FileData::getValuesAt(double x, double y, double *values)
{
    int i,j,k;
    bool valid = false;
    double theta = atan(y/x);
    double r = sqrt(x*x + y*y);

    // Find theta location
    if (theta < this->t_jph[0]) return;
    for (i=0; i<this->numTheta; ++i) {
        if (theta < this->t_jph[i+1]) {
            valid = true;
            break;
        }
    }
    if (!valid) return;

    // Find r location
    valid = false;
    if (r < this->r_iph[i][0]) return;
    for (j=0; j<this->Nr[i]-1; ++j) {
        if (r < this->r_iph[i][j+1]) {
            valid = true;
            break;
        }
    }
    if (!valid) return;

    // Updates the values
    for (k=0; k<this->numQ; ++k) {
        values[k] = this->cells[k][i][j+1];
    }
}

int FileData::getValuesAtTheta(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)
{
    int i,j;
    bool valid = false;
    double theta = atan(y/x);

    if (theta < this->t_jph[0]) return 0;
    for (i=0; i<this->numTheta; ++i) {
        if (theta < this->t_jph[i+1]) {
            valid = true;
            break;
        }
    }
    if (!valid) return 0;

    // Populate x & y data
    for (j=0; j<this->Nr[i]; ++j) {
        xdata->push_back(this->r_iph[i][j]);
        ydata->push_back(this->cells[n][i][j]);
    }

    return Nr[i];
}

int FileData::getValuesAtR(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)
{
    int i,j, count=0;
    bool valid;
    double r = sqrt(x*x + y*y);

    for (i=0; i<this->numTheta; ++i) {
        valid = false;
        if (r < this->r_iph[i][0]) continue;
        for (j=0; j<this->Nr[i]-1; ++j) {
            if (r < this->r_iph[i][j+1]) {
                valid = true;
                break;
            }
        }
        if (!valid) continue;

        count ++;
        xdata->push_back(this->t_jph[i]);
        ydata->push_back(this->cells[n][i][j]);
    }

    return count;
}

void FileData::readPatch(hid_t *h5dst, void *data, hid_t type, int *start, int *loc_size, int *glo_size)
{
    int dim = 2, d;

    hsize_t mdims[dim];
    hsize_t fdims[dim];

    hsize_t fstart[dim];
    hsize_t fstride[dim];
    hsize_t fcount[dim];
    hsize_t fblock[dim];

    for (d=0 ; d<dim ; ++d) {
        mdims[d] = loc_size[d];
        fdims[d] = glo_size[d];

        fstart[d]  = start[d];
        fstride[d] = 1;
        fcount[d]  = loc_size[d];
        fblock[d]  = 1;
    }

    hid_t mspace = H5Screate_simple(dim,mdims,NULL);
    hid_t fspace = H5Screate_simple(dim,fdims,NULL);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, fstart, fstride, fcount, fblock);

    H5Dread(*h5dst, type, mspace, fspace, H5P_DEFAULT, data);

    H5Sclose(mspace);
    H5Sclose(fspace);
}

// In a surprising twist, freeAll does not actually free
// this->minmax.
void FileData::freeAll()
{
    if ( !this->fileLoaded ) return;

    if (this->Nr != NULL) free(this->Nr);
    if (this->t_jph != NULL) free(this->t_jph);

    if (this->minmax != NULL) {
        if (this->minmax[0] != NULL) free(this->minmax[0]);
        if (this->minmax[1] != NULL) free(this->minmax[1]);
    }

    int i,j;
    for (i=0; i<this->numTheta; ++i) {
        if (this->r_iph != NULL && this->r_iph[i] != NULL) free(this->r_iph[i]);
        if (this->cells != NULL) {
            for (j=0; j<this->numQ; ++j) {
                if (this->cells[j][i] != NULL) free(this->cells[j][i]);
            }
        }
    }

    if (this->r_iph != NULL) free(this->r_iph);
    if (this->cells != NULL) {
        for (j=0; j<this->numQ; ++j) if (this->cells[j] != NULL) free(this->cells[j]);
        free(this->cells);
    }

}
