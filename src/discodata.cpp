#include <float.h>
#include "filedata.h"
#include "discodata.h"

DiscoData::DiscoData()
{
    this->numCells = 0;
    this->numR = 1;
    this->numZ = 1;
    this->time = 0.0;
    this->fileLoaded = false;

    this->minmax = (double **)malloc(sizeof(double *) * 2);

    // Initialize null pointers;
    this->cells = NULL;
    this->p_iph = NULL;
    this->r_jph = NULL;
    this->z_kph = NULL;
    this->Np = NULL;
    this->minmax[0] = NULL;
    this->minmax[1] = NULL;
    this->bounds[0] = 0.0;
    this->bounds[1] = 1.0;
    this->bounds[2] = 0.0;
    this->bounds[3] = 1.0;
}

DiscoData::~DiscoData()
{
    this->freeAll();
}

bool DiscoData::loadFromFile(const char *filename)
{
    // Useful variables
    int i, j, q, count;
    hsize_t dims[3];
    int start[2];
    int loc_size[2];
    int glo_size[2];
    int *pindex;

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

        // Get number of z plates
        hid_t h5_z = H5Dopen1(h5gridgroup, "z_kph");
        hid_t h5_z_spc = H5Dget_space(h5_z);
        H5Sget_simple_extent_dims(h5_z_spc, dims, NULL);
        this->numZ = dims[0] - 1; // TODO ... setting things here causes issues. ???
        H5Sclose(h5_z_spc);
        
        // Then we actually load the z_kph heights
        this->z_kph = (double *)malloc(dims[0] * sizeof(double));
        start[0] = 0; start[1] = 0;
        loc_size[0] = dims[0]; loc_size[1] = 1;
        glo_size[0] = dims[0]; glo_size[1] = 1;
        readPatch(&h5_z, this->z_kph, H5T_NATIVE_DOUBLE, start, loc_size, glo_size);
        H5Dclose(h5_z);
    
        // Get number of r annuli
        hid_t h5_r = H5Dopen1(h5gridgroup, "r_jph");
        hid_t h5_r_spc = H5Dget_space(h5_r);
        H5Sget_simple_extent_dims(h5_r_spc, dims, NULL);
        this->numR = dims[0] - 1; // TODO ... setting things here causes issues. ???
        H5Sclose(h5_r_spc);

        // Then we actually load the r_jph radii
        this->r_jph = (double *)malloc(dims[0] * sizeof(double));
        start[0] = 0; start[1] = 0;
        loc_size[0] = dims[0]; loc_size[1] = 1;
        glo_size[0] = dims[0]; glo_size[1] = 1;
        readPatch(&h5_r, this->r_jph, H5T_NATIVE_DOUBLE, start, loc_size, glo_size);
        H5Dclose(h5_r);

        //We only load the equatorial slice!
        int k = this->numZ / 2;

        // ... and we get the number of phi cells per annulus
        this->Np = (int *)malloc(this->numR * sizeof(int));
        hid_t h5_Np = H5Dopen1(h5gridgroup, "Np");
        start[0] = 0; start[1] = k;
        loc_size[0] = this->numR; loc_size[1] = 1;
        glo_size[0] = this->numR; glo_size[1] = 1;
        readPatch(&h5_Np, this->Np, H5T_NATIVE_INT, start, loc_size, glo_size);
        H5Dclose(h5_Np);

        // ... and finally, we load the indices for the cells, into a temp array.
        pindex = (int *)malloc(this->numR * sizeof(int));
        hid_t h5_pindex = H5Dopen1(h5gridgroup, "Index");
        readPatch(&h5_pindex, pindex, H5T_NATIVE_INT, start, loc_size, glo_size);
        H5Dclose(h5_pindex);

    // Then we load the information in the cells (geometry + values)
    hid_t h5_cells = H5Dopen1(h5datagroup, "Cells");
    hid_t h5_cells_spc = H5Dget_space(h5_cells);
    H5Sget_simple_extent_dims(h5_cells_spc, dims, NULL);

    this->numCells = 0;
    for(j=0; j<this->numR; j++)
        this->numCells += this->Np[j];
    this->numQ = dims[1] - 1;
    H5Sclose(h5_cells_spc);

    // Set up minmax
    this->minmax[0] = (double *)malloc(sizeof(double) * this->numQ);
    this->minmax[1] = (double *)malloc(sizeof(double) * this->numQ);

    this->p_iph = (double **)malloc(sizeof(double *) * this->numR);
    this->cells = (double **)malloc(sizeof(double *) * this->numQ);
    for (i=0; i<this->numQ; ++i) {
        cells[i] = (double *)malloc(sizeof(double) * this->numCells);
    }

    loc_size[1] = this->numQ+1;
    glo_size[0] = dims[0]; glo_size[1] = this->numQ + 1;
    start[1] = 0;
    
    for (q=0; q<this->numQ; q++) {
        this->minmax[0][j] = DBL_MAX;
        this->minmax[1][j] = -DBL_MAX;
    }
    
    count = 0;
    for (j=0; j<this->numR; j++)
    {
        this->p_iph[j] = (double *)malloc(sizeof(double) * Np[j]);
        loc_size[0] = this->Np[j];
        start[0] = pindex[j];
        double buffer[this->Np[j] * (this->numQ+1)];
        readPatch(&h5_cells, buffer, H5T_NATIVE_DOUBLE, start, loc_size, glo_size);

        // Load from the buffer into the proper class variables
        for (i=0; i<Np[j]; i++)
        {
            this->p_iph[j][i] = buffer[i * (this->numQ + 1) + this->numQ];
            for (q = 0; q < this->numQ; q++)
            {
                this->cells[q][count+i] = buffer[i * (this->numQ + 1) + q];
                if (this->minmax[0][q] > this->cells[k][count+j]) 
                    this->minmax[0][q] = this->cells[k][count+j];
                if (this->minmax[1][q] < this->cells[k][count+j]) 
                    this->minmax[1][q] = this->cells[k][count+j];
            }
        }
        count += Np[j];
    }

    this->bounds[0] = -this->r_jph[this->numR];
    this->bounds[1] =  this->r_jph[this->numR];
    this->bounds[2] = -this->r_jph[this->numR];
    this->bounds[3] =  this->r_jph[this->numR];

    // Housekeeping!
    H5Dclose(h5_cells);
    H5Gclose(h5gridgroup);
    H5Gclose(h5datagroup);
    H5Fclose(h5file);

    free(pindex);

    this->fileLoaded = true;

    #if FILEDATA_VERBOSE
    fprintf(stderr, "Done loading data.\n");
    #endif

    return true;
}

// This function sets values[i] to the i'th q of the cell
// located at (x,y). If such a cell does not exist, it
// simply returns.
void DiscoData::getValuesAt(double x, double y, double *values)
{
    int i,j,q;
    bool valid = false;
    bool before = false;
    double phi = atan2(y,x);
    double r = sqrt(x*x + y*y);

    // Find r location
    int count = 0;
    if(r < this->r_jph[0])
        return;
    for(j = 0; j < this->numR; j++)
    {
        if (r < this->r_jph[j+1]) 
        {
            valid = true;
            break;
        }
        count += this->Np[j];
    }
    if (!valid) 
        return;

    // Find phi location
    valid = false;

    for(i = 0; i < this->Np[j]; i++)
    {
        double diff = this->p_iph[j][i] - phi;
        while(diff < -M_PI)
            diff += 2*M_PI;
        while(diff > M_PI)
            diff -= 2*M_PI;

        if (diff < 0.0) 
            before = true;
        if(before && diff > 0.0)
        {
            valid = true;
            break;
        }
    }
    if(!valid)
    {
        i = 0;
        valid = true;
    }

    // Updates the values
    for (q=0; q<this->numQ; q++)
        values[q] = this->cells[q][count+i];
}

int DiscoData::getValuesAtTheta(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)
{
    int i,j, size=0, count=0;
    bool valid;
    double phi = atan2(y,x);

    for(j = 0; j < this->numR; j++)
    {
        valid = false;

        for(i = 0; i < this->Np[j]; i++)
        {
            double diff = this->p_iph[j][i] - phi;
            while(diff < -M_PI)
                diff += 2*M_PI;
            while(diff > M_PI)
                diff -= 2*M_PI;

            if (diff < 0.0) 
                valid = true;
            if(valid && diff > 0.0)
                break;
        }
        if(!valid)
            i = 0;

        xdata->push_back(0.5*(this->r_jph[j]+this->r_jph[j+1]));
        ydata->push_back(this->cells[n][count+i]);

        count += this->Np[j];
    }

    size = this->numR;

    return size;
}

int DiscoData::getValuesAtR(double x, double y, int n, std::vector<double> *xdata, std::vector<double> *ydata)
{
    int i,j,count=0;
    bool valid = false;
    double r = sqrt(x*x+y*y);

    if (r < this->r_jph[0])
        return 0;
    for (j=0; j<this->numR; j++)
    {
        if (r < this->r_jph[j+1]) {
            valid = true;
            break;
        }
        count += this->Np[j];
    }
    if (!valid) return 0;

    // Populate x & y data
    for (i=0; i<this->Np[j]; i++)
    {
        xdata->push_back(0.5*(this->p_iph[j][i]+this->p_iph[j][i+1]));
        ydata->push_back(this->cells[n][count+i]);
    }

    return Np[j];
}

void DiscoData::genGridData(double **gp, int *ngp, int **ci, int *nci, int **gpi, int *ngpi)
{
    int i,j,c,l,v;

    #if FILEDATA_VERBOSE
    fprintf(stderr, "Generating grid data...\n");
    #endif

    *nci = 4*this->numCells;
    *ngp = 4*this->numCells;
    *ngpi = 4*this->numCells;

    *gp = (double *) malloc(2 * (*ngp) * sizeof(double));
    *ci = (int *) malloc((*nci) * sizeof(int));
    *gpi = (int *) malloc((*ngpi) * sizeof(int));

    double p0,p1,dp, r0,r1,dr;

    l=0;
    v=0;
    
    for(j=0; j<this->numR; j++)
    {
        r0 = this->r_jph[j]; 
        r1 = this->r_jph[j+1];
        dr = (r1 - r0) * 0.01;
        r0 -= dr; r1 += dr;
        
        for (i=0; i<this->Np[j]; i++)
        {

            // DEBUG

            p0 = this->p_iph[j][i];
            if(i+1 < this->Np[j])
                p1 = this->p_iph[j][i+1];
            else
                p1 = this->p_iph[j][0];
            dp = (p1 - p0);
            while(dp < 0.0)
                dp += 2*M_PI;
            while(dp > 2*M_PI)
                dp -= 2*M_PI;
            dp *= 0.05;
            p0 -= dp; p1 += dp;

            // Bottom right
            (*gp)[2*v] = r0*cos(p0);
            (*gp)[2*v+1] = r0*sin(p0);
            (*ci)[v] = v;
            v++;

            // Bottom left
            (*gp)[2*v] = r0*cos(p1);
            (*gp)[2*v+1] = r0*sin(p1);
            (*ci)[v] = v;
            v++;

            // Top left
            (*gp)[2*v] = r1*cos(p1);
            (*gp)[2*v+1] = r1*sin(p1);
            (*ci)[v] = v;
            v++;

            // Top right
            (*gp)[2*v] = r1*cos(p0);
            (*gp)[2*v+1] = r1*sin(p0);
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
    
    #if FILEDATA_VERBOSE
    fprintf(stderr, "Done generating grid data.\n");
    #endif
}

//Free's only DiscoData specific elements. FileData destructor takes care
//of the rest.
void DiscoData::freeAll()
{
    if ( !this->fileLoaded ) 
        return;

    if (this->Np != NULL) 
        free(this->Np);
    if (this->r_jph != NULL) 
        free(this->r_jph);
    if (this->z_kph != NULL) 
        free(this->z_kph);

    int j;
    for (j=0; j<this->numR; ++j) {
        if (this->p_iph != NULL && this->p_iph[j] != NULL) 
            free(this->p_iph[j]);
    }
    if (this->p_iph != NULL)
        free(this->p_iph);
}

