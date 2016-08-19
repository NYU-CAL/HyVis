#include "filedata.h"

FileData::FileData()
{
    this->numCells = 0;
    this->time = 0.0;
    this->fileLoaded = false;

    this->minmax = (double **)malloc(sizeof(double *) * 2);

    // Initialize null pointers;
    this->cells = NULL;
    this->minmax[0] = NULL;
    this->minmax[1] = NULL;
    this->bounds[0] = 0.0;
    this->bounds[1] = 1.0;
    this->bounds[2] = 0.0;
    this->bounds[3] = 1.0;
}

FileData::~FileData()
{
    free(this->minmax);
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

int FileData::get_nq()
{
    if (this->fileLoaded) return this->numQ;
    return 0;
}

double **FileData::get_cells()
{
    if (this->fileLoaded) return this->cells;
    return NULL;
}

double *FileData::get_bounds()
{
    if (this->fileLoaded) return this->bounds;
    return NULL;
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
    int j;

    if ( !this->fileLoaded ) return;

    if (this->minmax != NULL) {
        if (this->minmax[0] != NULL) free(this->minmax[0]);
        if (this->minmax[1] != NULL) free(this->minmax[1]);
    }

    if (this->cells != NULL) {
        for (j=0; j<this->numQ; ++j) 
            if (this->cells[j] != NULL) 
                free(this->cells[j]);
        free(this->cells);
    }
}

