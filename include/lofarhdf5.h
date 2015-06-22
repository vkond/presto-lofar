#ifndef _LOFARHDF5_H
#define _LOFARHDF5_H

using namespace std;
// we need this undef, because there is #define READWRITE in fitsio.h
// used by psrfits, and they are in conflict, because #define is global
#undef READWRITE
#include "dal/lofar/BF_File.h"
using namespace dal;

#ifdef  __cplusplus
extern "C" {
#endif

// struct to read input data in HDF5 format 
typedef struct {
  BF_File* fd; // file descriptor
  BF_StokesDataset* bf_stokes;
  long long current_sample;
} h5file;

#ifdef  __cplusplus
}
#endif

/* lofarhdf5.cxx */
long long offset_to_LOFARHDF5_spectra(long long specnum, struct spectra_info *s);
int get_LOFARHDF5_rawblock(float *fdata, struct spectra_info *s, int *padding);

#endif /* _LOFARHDF5_H */
