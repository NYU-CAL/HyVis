#include "filedata.h"
#include "jetdata.h"

JetData::JetData()
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
    this->bounds[0] = 0.0;
    this->bounds[1] = 1.0;
    this->bounds[2] = 0.0;
    this->bounds[3] = 1.0;
}

JetData::~JetData()
{
    this->freeAll();
}

bool JetData::loadFromFile(const char *filename)
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
    this->cells = (double **)malloc(sizeof(double *) * this->numQ);
    for (i=0; i<this->numQ; ++i) {
        cells[i] = (double *)malloc(sizeof(double) * this->numCells);
    }

    loc_size[1] = this->numQ+1;
    glo_size[0] = this->numCells; glo_size[1] = this->numQ + 1;
    start[0] = - this->Nr[0];
    
    int count = 0;
    for (i=0; i<this->numTheta; ++i) {
        this->r_iph[i] = (double *)malloc(sizeof(double) * (Nr[i]+1));
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
        this->r_iph[i][0] = 0.0;
        for (j=0; j<Nr[i]; ++j) {
            this->r_iph[i][j+1] = buffer[j * (this->numQ + 1) + this->numQ];
            for (int k=0; k<this->numQ; ++k) {
                this->cells[k][count+j] = buffer[j * (this->numQ + 1) + k];
                if (this->minmax[0][k] > this->cells[k][count+j]) this->minmax[0][k] = this->cells[k][count+j];
                if (this->minmax[1][k] < this->cells[k][count+j]) this->minmax[1][k] = this->cells[k][count+j];
            }
        }
        count += Nr[i];
    }

    double th0 = this->t_jph[0];
    double th1 = this->t_jph[this->numTheta];
    double r0 = this->r_iph[0][0];
    double r1 = this->r_iph[0][this->Nr[0]];
    double xmin = fmin(r0*cos(th1), r1*cos(th1));
    double xmax = r1*cos(th0);
    double ymin = fmin(r0*sin(th0), r0*sin(th1));
    double ymax = fmax(r1*sin(th0), r1*sin(th1));
    this->bounds[0] = xmin;
    this->bounds[1] = xmax;
    this->bounds[2] = ymin;
    this->bounds[3] = ymax;

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

// This function sets values[i] to the i'th q of the cell
// located at (x,y). If such a cell does not exist, it
// simply returns.
void JetData::getValuesAt(double x, double y, double *values)
{
    int i,j,k;
    bool valid = false;
    double theta = atan(y/x);
    double r = sqrt(x*x + y*y);

    // Find theta location
    int count = 0;
    if (theta < this->t_jph[0]) return;
    for (i=0; i<this->numTheta; ++i) {
        if (theta < this->t_jph[i+1]) {
            valid = true;
            break;
        }
        count += this->Nr[i];
    }
    if (!valid) return;

    // Find r location
    valid = false;
    if (r < this->r_iph[i][0]) return;
    for (j=0; j<this->Nr[i]; ++j) {
        if (r < this->r_iph[i][j+1]) {
            valid = true;
            break;
        }
    }
    if (!valid) return;

    // Updates the values
    for (k=0; k<this->numQ; ++k) {
        values[k] = this->cells[k][count+j];
    }
}

int JetData::getValuesAtTheta(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)
{
    int i,j,count=0;
    bool valid = false;
    double theta = atan(y/x);

    if (theta < this->t_jph[0]) return 0;
    for (i=0; i<this->numTheta; ++i) {
        if (theta < this->t_jph[i+1]) {
            valid = true;
            break;
        }
        count += this->Nr[i];
    }
    if (!valid) return 0;

    // Populate x & y data
    for (j=0; j<this->Nr[i]; ++j) {
        xdata->push_back(0.5*(this->r_iph[i][j]+this->r_iph[i][j+1]));
        ydata->push_back(this->cells[n][count+j]);
    }

    return Nr[i];
}

int JetData::getValuesAtR(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)
{
    int i,j, size=0, count=0;
    bool valid;
    double r = sqrt(x*x + y*y);

    for (i=0; i<this->numTheta; ++i) {
        valid = false;
        if (r < this->r_iph[i][0]) continue;
        for (j=0; j<this->Nr[i]; ++j) {
            if (r < this->r_iph[i][j+1]) {
                valid = true;
                break;
            }
        }
        if (!valid) continue;

        size ++;
        xdata->push_back(0.5*(this->t_jph[i]+this->t_jph[i+1]));
        ydata->push_back(this->cells[n][count+j]);

        count += this->Nr[i];
    }

    return size;
}

void JetData::genGridData(double **gp, int *ngp, int **ci, int *nci, int **gpi, int *ngpi)
{
    int i,j,c,l,v;

    *nci = 4*this->numCells;
    *ngp = 4*this->numCells;
    *ngpi = 4*this->numCells;

    *gp = (double *) malloc(2 * (*ngp) * sizeof(double));
    *ci = (int *) malloc((*nci) * sizeof(int));
    *gpi = (int *) malloc((*ngpi) * sizeof(int));

    double t0,t1,dt, r0,r1,dr;

    l=0;
    v=0;

    if (this->numTheta > 0) {
        for (i=0; i<this->numTheta; ++i) {
            t0 = this->t_jph[i]; 
            t1 = this->t_jph[i+1];
            dt = (t1 - t0) * 0.01;
            t0 -= dt; t1 += dt;

            for (j=0; j<this->Nr[i]; ++j) {

                // DEBUG

                r0 = this->r_iph[i][j];
                r1 = this->r_iph[i][j+1];
                dr = (r1 - r0) * 0.05;
                r0 -= dr; r1 += dr;

                //r0 = r_iph[0][j-1]; r1 = r_iph[0][j];
                

                // Bottom right
                (*gp)[2*v] = r0*cos(t0);
                (*gp)[2*v+1] = r0*sin(t0);
                (*ci)[v] = v;
                v++;

                // Bottom left
                (*gp)[2*v] = r0*cos(t1);
                (*gp)[2*v+1] = r0*sin(t1);
                (*ci)[v] = v;
                v++;

                // Top left
                (*gp)[2*v] = r1*cos(t1);
                (*gp)[2*v+1] = r1*sin(t1);
                (*ci)[v] = v;
                v++;

                // Top right
                (*gp)[2*v] = r1*cos(t0);
                (*gp)[2*v+1] = r1*sin(t0);
                (*ci)[v] = v;
                v++;

                // BR -> TR
                (*gpi)[l] = v-4; l++;
                (*gpi)[l] = v-1; l++;

                // BR -> BL
                (*gpi)[l] = v-4; l++;
                (*gpi)[l] = v-3; l++;
            }
        }
    } else {
        for (j=0; j<this->Nr[0]; ++j) {
            r0 = r_iph[0][j]; r1 = r_iph[0][j+1];
            dr = (r1 - r0) * 0.05;
            r0 -= dr; r1 += dr;

            (*ci)[c] = c; c++;

            // Bottom right
            (*gp)[2*v] = r0;
            (*gp)[2*v+1] = 1.0;
            v++;

            // Bottom left
            (*gp)[2*v] = r0;
            (*gp)[2*v+1] = 1.0;
            v++;

            // Top left
            (*gp)[2*v] = r1;
            (*gp)[2*v+1] = 1.0;
            v++;

            // Top right
            (*gp)[2*v] = r1;
            (*gp)[2*v+1] = 0.0;
            v++;

            // BR -> TR
            (*gpi)[l] = v-4; l++;
            (*gpi)[l] = v-1; l++;

            // BR -> BL
            (*gpi)[l] = v-4; l++;
            (*gpi)[l] = v-3; l++;
        }
    }

}

//Free's only JetData specific elements. FileData destructor takes care
//of the rest.
void JetData::freeAll()
{
    if ( !this->fileLoaded ) 
        return;

    if (this->Nr != NULL) 
        free(this->Nr);
    if (this->t_jph != NULL) 
        free(this->t_jph);

    int i,j;
    for (i=0; i<this->numTheta; ++i) {
        if (this->r_iph != NULL && this->r_iph[i] != NULL) 
            free(this->r_iph[i]);
    }

    if (this->r_iph != NULL) free(this->r_iph);
}

