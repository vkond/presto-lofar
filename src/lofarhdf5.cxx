extern "C" {
    #include "presto.h"
}
#include "lofarhdf5.h"

#include "backend_common.h"

using namespace std;
// we need this undef, because there is #define READWRITE in fitsio.h
// used by psrfits, and they are in conflict, because #define is global
#undef READWRITE
#include "dal/lofar/BF_File.h"
using namespace dal;

#include <iomanip>  // for setprecision use
#include <iostream>

static h5file *h5s; // array of h5file structs

static int currentfile = 0, currentblock = 0;
static int numbuffered = 0, numpadded = 0;

extern "C" void add_padding(float *fdata, float *padding, int numchan, int numtopad);

long long offset_to_LOFARHDF5_spectra(long long specnum, struct spectra_info *s);
int get_LOFARHDF5_rawblock(float *fdata, struct spectra_info *s, int *padding);



/* Reads LOFAR HDF5 input data */
/* Only supports for now Stokes I single beam data */
/* All subbands should be in one file */
void read_LOFARHDF5_files(struct spectra_info *s) {

  // e.g. for debugging purposes
  int verbose = 0;

  s->datatype = LOFARHDF5;
  if (s->num_files > 1) {
      cerr << "WARNING: Only 1st input file will be processed.\n No support of multiple LOFAR HDF5 input files yet" << endl;
      s->num_files = 1;
  }
  s->header_offset = gen_ivect(s->num_files);
  s->start_spec = (long long *)malloc(sizeof(long long) * s->num_files);
  s->num_spec = (long long *)malloc(sizeof(long long) * s->num_files);
  s->num_pad = (long long *)malloc(sizeof(long long) * s->num_files);
  s->start_MJD = (long double *)malloc(sizeof(long double) * s->num_files);

  BF_File *fd = new BF_File (s->filenames[0]);
  cerr << "Reading LOFAR HDF5 file..." << endl << endl;

  // getting PI
  Attribute<std::string> observer = fd->projectPI();
  if (observer.exists()) {
      strcpy(s->observer, observer.get().c_str());
      if (verbose) cerr << "PI=" << s->observer << endl;
  } else {
    // getting project contact
    Attribute<std::string> projectContact = fd->projectContact();
    if (projectContact.exists()) {
        strcpy(s->observer, projectContact.get().c_str());
        if (verbose) cerr << "PROJECT_CONTACT=" << s->observer << endl;
    }
  }

  // getting project id
  Attribute<std::string> obsid = fd->observationID();
  if (obsid.exists()) {
    strcpy(s->project_id, obsid.get().c_str());
    if (verbose) cerr << "ObsID=" << s->project_id << endl;
  }

  // getting start UTC time
  Attribute<std::string> startUTC = fd->observationStartUTC();
  if (startUTC.exists()) {
    strcpy(s->date_obs, startUTC.get().c_str());
    if (verbose) cerr << "Start UTC=" << s->date_obs << endl;
  }

  // getting targets
  Attribute< std::vector<std::string> > BFtargets = fd->targets();
  if (verbose) {
   if (BFtargets.exists()) {
      std::vector<std::string> t = BFtargets.get();
      std::vector<std::string>::size_type i=0;
      for(std::vector<std::string>::iterator it = t.begin(); it != t.end(); ++it, i++)
        cerr << "TARGET" << i << "=" << *it << endl;
   } else cerr << "TARGET does not exist" << endl;
  }

  // getting frequency center
  Attribute<double> freq = fd->observationFrequencyCenter();
  if (!freq.exists()) {
   cerr << "observationFrequencyCenter not defined" << endl;
   exit (1);
  } else { 
    s->fctr = freq.get();
    if (verbose) cerr << "observation frequency=" << setprecision(20) << s->fctr << " MHz" << endl;
  }

  // getting number of SAPs
  Attribute<unsigned> nsap = fd->observationNofSubArrayPointings();
  if (!nsap.exists()) {
   cerr << "observationNofSubArrayPointings not defined" << endl;
   exit (1);
  } else if (verbose) cerr << "number of SAPs=" << nsap.get() << endl;

  // getting  the instance of SAP
  // checking all SAPs if they exist to pick the right one (if there will be two SAPs in one file, only
  // the first one will be picked up)
  unsigned sap_index;
  for (sap_index=0; sap_index<nsap.get(); sap_index++) {
    if (fd->subArrayPointing(sap_index).exists()) break;
  }
  BF_SubArrayPointing sap = fd->subArrayPointing(sap_index);

  Attribute<unsigned> nbeam = sap.observationNofBeams();
  if (!nbeam.exists()) {
   cerr << "sap.observationNofBeams not defined" << endl;
   exit (1);
  } else if (verbose) cerr << "number of beams=" << nbeam.get() << endl;

  // getting the instance of first TA beam in the SAP
  // checking all TABs in the SAP if they exist in the file until the first one that exists is found
  unsigned tab_index;
  for (tab_index=0; tab_index<nbeam.get(); tab_index++) {
    if (sap.beam(tab_index).exists()) break;
  }
  BF_BeamGroup beam = sap.beam(tab_index);

  // getting the center frequency of the beam
  Attribute<double> freq2 = beam.beamFrequencyCenter();
  if (!freq2.exists()) {
   cerr << "beam.beamFrequencyCenter not defined" << endl;
   exit (1);
  } else if (verbose) cerr << "beam frequency=" << setprecision(20) << freq2.get() << endl;

  // getting the subband width
  Attribute<double> bw2 = beam.subbandWidth();
  if (!bw2.exists()) {
   cerr << "beam.subbandWidth not defined" << endl;
   exit (1);
  } else if (verbose) cerr << "sap subbandwidth=" << setprecision(20) << bw2.get() << endl;

  // getting number of channels per subband
  Attribute<unsigned> nchan = beam.channelsPerSubband();
  if (!nchan.exists()) {
   cerr << "beam.channelsPerSubband not defined" << endl;
   exit (1);
  } else if (verbose) cerr << "number of channels/sub=" << nchan.get() << endl;

  // getting the pointer for the Stokes class
  BF_StokesDataset *bf_stokes = NULL;

  for (unsigned i=0; i<4; i++) {
    BF_StokesDataset tmp = beam.stokes(i);
    if (tmp.exists()) {
      bf_stokes = new BF_StokesDataset (beam.stokes(i));
      break;
    }
  }
  if (bf_stokes == NULL) {
    cerr << "No BF_StokesDataset present" << endl;
    exit (1);
  }

  // getting the Stokes component
  Attribute<std::string> stokesC = bf_stokes->stokesComponent();
  if (verbose) if (stokesC.exists()) cerr << "stokes component=" << stokesC.get() << endl;

  // getting the number of subbands
  Attribute<unsigned> nsub = bf_stokes->nofSubbands();
  if (nsub.exists()) {
   if (verbose) cerr << "nsub=" << nsub.get() << endl;
  } else { if (verbose) cerr << "stokes nofSubbands not defined" << endl; }

  // getting the number of channels for each subband
  Attribute< std::vector<unsigned> > nofchan = bf_stokes->nofChannels();
  if (verbose) if (nchan.exists()) {
                 std::vector<unsigned> nchan = nofchan.get();
                 cerr << "stokes nofChannels size=" << nchan.size() << endl;
               } else cerr << "stokes nofChannels not defined" << endl;

  // getting the rank of the dataset
  size_t ndim= bf_stokes->ndims();
  if (verbose) cerr << "stokes ndim=" << ndim << endl;

  if (verbose) {
   std::vector<std::string> files = bf_stokes->externalFiles();
   for (unsigned i=0; i<files.size(); i++)
     cerr << "files[" << i << "]=" << files[i] << endl;
  }

  // getting telescope
  Attribute<std::string> telescope = fd->telescope();
  if (telescope.exists()) {
   strcpy(s->telescope, telescope.get().c_str());
   if (verbose) cerr << "telescope=" << s->telescope << endl;
  }

  // setting machine
  // For now assuming it is LOFAR's COBALT
  strcpy(s->backend, "COBALT");
  // setting frontend
  Attribute<std::string> filter = fd->filterSelection();
  if (filter.exists()) {
    strcpy(s->frontend, filter.get().c_str());
    if (verbose) cerr << "FRONTEND=" << s->frontend << endl;
  }

  // getting the vector of targets
  Attribute< std::vector<std::string> > targets = beam.targets();
  if (targets.exists()) {
    std::vector<std::string> t = targets.get();
    if (t.size() != 0) {
     strcpy(s->source, t.front().c_str());
     if (verbose) cerr << "target = " << t.front() << endl;
    } else { if (verbose) cerr << "targets vector is empty" << endl; }
  } else { if (verbose) cerr << "beam target does not exist" << endl; }

  // getting number of samples
  Attribute<unsigned> nsamp = bf_stokes->nofSamples();
  if (nsamp.exists()) s->N = nsamp.get();

  // are data in Complex Voltage format?
  Attribute<bool> volts = beam.complexVoltage();
  if (volts.exists() && volts.get() == 1) {
    cerr << "Can't process complex-voltage data, ndim = " << ndim << endl;
    exit (1);
  }
 
  // check for which coordinate is Spectral
  unsigned spectral_dim = 1;

  // getting instance of Coordinates container
  CoordinatesGroup coord = beam.coordinates();
  if (coord.exists()) {
    Attribute< std::vector<std::string> > types = coord.coordinateTypes();
    if (types.exists()) {
      std::vector<std::string> t = types.get();
      for (unsigned i=0; i<t.size(); i++) {
	if (t[i] == "Spectral") {
	  spectral_dim = i;
	  break;
	}
      }
    }
  }

  std::vector<ssize_t> dims = bf_stokes->dims();
  s->orig_num_chan = dims[spectral_dim];
  s->num_channels = s->orig_num_chan;
  if (verbose) cerr << "Total number of channels=" << s->orig_num_chan << endl;
  
  // getting number of Stokes components in one file
  //Attribute<unsigned> npol = beam.nofStokes();
  // getting number of Stokes components in the observation
  Attribute<unsigned> npol = beam.observationNofStokes();
  unsigned stokes_npol = 1;

  if (npol.exists()) stokes_npol = npol.get();

  if (stokes_npol == 1) {
    s->summed_polns = 1;
    s->num_polns = 1;
  } else {
    cerr << "Can't process more than one IFs" << endl;
    exit (1);
  }

  s->bits_per_sample = 32;
  s->num_beams = 1; // For now assuming there is only one beam
  s->beamnum = 0;

  // getting split Frequency center of the beam
  Attribute<double> cfreq = beam.beamFrequencyCenter();
  if (!cfreq.exists()) {
   cerr << "beamFrequencyCenter not defined" << endl;
   exit (1);
  } else { if (verbose) cerr << "beamFrequencyCenter=" << setprecision(20) << cfreq.get() << endl;
  }

  // getting the start MJD
  Attribute<double> mjd = fd->observationStartMJD();
  if (mjd.exists()) s->start_MJD[0] = mjd.get();
  if (verbose) cerr << "MJD=" << setprecision(20) << s->start_MJD[0] << endl;

  // getting the clock rate
  Attribute<double> cRate = fd->clockFrequency();
  if (verbose) {
    if (cRate.exists()) cerr << "clockRate=" << setprecision(20) << cRate.get() << endl;
    else cerr << "clockRate undefined" << endl;
  }

  // getting the sampling rate
  Attribute<double> sRate = beam.samplingRate();
  if (verbose) {
    if (sRate.exists()) cerr << "samplingRate=" << setprecision(20) << sRate.get() << endl;
    else cerr << "samplingRate undefined" << endl;
  }

  // getting the sampling time
  Attribute<double> sTime = beam.samplingTime();
  if (sTime.exists()) {
    s->dt = sTime.get();  
    if (verbose) cerr << "samplingTime=" << setprecision(20) << s->dt << " s"<< endl;
  } else if (verbose) cerr << "samplingTime undefined" << endl;

  // getting the channel width
  Attribute<double> rate = beam.channelWidth();
  if (!rate.exists()) {
   cerr << "beam.channelWidth not defined" << endl;
   exit (1);
  } else { s->orig_df = rate.get() * 1.0e-6;
           s->df = s->orig_df;
           if (verbose) cerr << "channel Width=" << setprecision(20) << s->orig_df << " Hz" << endl;
         }

  // getting the subband width
  Attribute<double> subwidth = beam.subbandWidth();
  if (verbose) if (subwidth.exists())
                 cerr << "subband Width=" << setprecision(20) << subwidth.get() << " Hz" << endl;
               else cerr << "subband Width undefined" << endl;

  // setting the bandwidth (in MHz) and center freqs of lowest and highest channels
  s->BW = s->orig_df * s->orig_num_chan;
  s->lo_freq = s->fctr - s->BW/2. + fabs(s->orig_df)/2.;
  s->hi_freq = s->fctr + s->BW/2. - fabs(s->orig_df)/2.;

  // getting the azimuth and zenith angle of the beam (in degrees)
  Attribute<std::vector <double> > azdeg = sap.pointAzimuth();
  if (azdeg.exists()) {
    std::vector<double> tmp = azdeg.get();
    std::vector<double>::iterator it = tmp.begin();
    s->azimuth = *it;
    if (verbose) cerr << "AZ=" << setprecision(20) << s->azimuth << " deg" << endl;
  } else if (verbose) cerr << "AZ undefined" << endl;

  Attribute<std::vector <double> > altdeg = sap.pointAltitude();
  if (altdeg.exists()) {
    std::vector<double> tmp = altdeg.get();
    std::vector<double>::iterator it = tmp.begin();
    s->zenith_ang = 90. - *it;
    if (verbose) cerr << "ZA=" << setprecision(20) << s->zenith_ang << " deg" << endl;
  } else if (verbose) cerr << "ZA undefined" << endl;

  // getting the RA and DEC of the beam (in degrees)
  Attribute<double> radeg = beam.pointRA();
  if (radeg.exists()) {
    s->ra2000 = radeg.get();
    if (verbose) cerr << "RA=" << setprecision(20) << s->ra2000 << " deg" << endl;
    int ra_h = (int)(s->ra2000/15.);
    int ra_m = (int)((s->ra2000/15. - ra_h)*60.);
    double ra_s = (s->ra2000/15. - ra_h - ra_m/60.)*3600.;
    sprintf(s->ra_str, "%02d:%02d:%s%lf", ra_h, ra_m, ra_s < 10 ? "0" : "", ra_s);
    if (verbose) cerr << "RA=" << s->ra_str << endl;
  } else if (verbose) cerr << "RA undefined" << endl;

  Attribute<double> decdeg = beam.pointDEC();
  if (decdeg.exists()) {
    s->dec2000 = decdeg.get();
    if (verbose) cerr << "DEC=" << setprecision(20) << s->dec2000 << " deg" << endl;
    int dec_d = (int)(fabs(s->dec2000));
    int dec_m = (int)((fabs(s->dec2000) - dec_d)*60.);
    double dec_s = (fabs(s->dec2000) - dec_d - dec_m/60.)*3600.;
    int sign = (int)(s->dec2000);
    if (sign < 0) dec_d = -dec_d;
    sprintf(s->dec_str, "%02d:%02d:%s%lf", dec_d, dec_m, dec_s < 10 ? "0" : "", dec_s);
    if (verbose) cerr << "DEC=" << s->dec_str << endl;
  } else if (verbose) cerr << "DEC undefined" << endl;

  // total time
  s->T = s->N * s->dt;  
  
  s->samples_per_spectra = s->num_polns * s->num_channels;
  s->bytes_per_spectra = s->bits_per_sample * s->samples_per_spectra / 8;
  s->spectra_per_subint = 480;  // same as for Sigproc filterbank data
  s->bytes_per_subint = s->bytes_per_spectra * s->spectra_per_subint;
  s->samples_per_subint = s->spectra_per_subint * s->samples_per_spectra;
  s->min_spect_per_read = 1;  // Can read a single spectra at a time
  s->apply_flipband = 0;

  s->time_per_subint = s->spectra_per_subint * s->dt;
  s->start_spec[0] = 0L;
  s->num_spec[0] = s->N;
  s->header_offset[0] = 0;
  s->get_rawblock = &get_LOFARHDF5_rawblock;
  s->offset_to_spectra = &offset_to_LOFARHDF5_spectra;

  s->h5files = (long long **)malloc(sizeof(long long*) * s->num_files);
  h5s = new h5file[s->num_files]; 
  h5s[0].fd = fd;
  h5s[0].bf_stokes = bf_stokes;
  h5s[0].current_sample = 0;
  s->h5files[0] = (long long*)(&h5s[0]);

  // allocate the raw data buffers
  s->padvals = gen_fvect(s->num_channels);
  for (int ii = 0 ; ii < s->num_channels ; ii++) s->padvals[ii] = 0.0;
}


// properly closing opened HDF5 file
void close_LOFARHDF5_file (long long *ptr) {
    h5file *h5 = *ptr; // assign address to h5file object
    if (h5s != NULL) delete (h5s);
    /*
    if (h5->bf_stokes != NULL) delete(h5->bf_stokes);
    if (h5->fd != NULL) {
        h5->fd->close();
        delete(h5->fd);
    }
    */
}


long long offset_to_LOFARHDF5_spectra(long long specnum, struct spectra_info *s)
// This routine offsets into the filterbank files to the spectra
// 'specnum'.  It returns the current spectra number.
{
    long long offset_spec;
    h5file *h5 = (h5file *)(s->h5files[0]);
    
    if (specnum > s->N) {
        fprintf(stderr, "Error:  offset spectra %lld is > total spectra %lld\n\n", 
               specnum, s->N);
        exit(1);
    }

    // For LOFAR data we have only one input file
    currentfile = 0;

    // Are we in a padding zone?
    if (specnum > (s->start_spec[currentfile] + s->num_spec[currentfile])) {
        // Seek to the end of the file
        // no need to do it actually for HDF5 files?...
        h5->current_sample = s->N;
        numbuffered = 0;
        numpadded = specnum - (s->start_spec[currentfile] + s->num_spec[currentfile]);
        return specnum;
    }

    // Otherwise, seek to the spectra
    // Skip the header first
    offset_spec = specnum - s->start_spec[currentfile];
    h5->current_sample = offset_spec;
    numbuffered = 0;
    numpadded = 0;
    return specnum;
}


int get_LOFARHDF5_rawblock(float *fdata, struct spectra_info *s, int *padding)
// This routine reads a single block (i.e subint) from the input files
// which contain raw data in LOFAR HDF5 format.  If padding is
// returned as 1, then padding was added and statistics should not be
// calculated.  Return 1 on success.
{
    int numread, numtopad = 0, numtoread;
    float *fdataptr = fdata;
    h5file *h5 = (h5file *)s->h5files[0];
    
    // If we only made a partial read last time, adjust our pointer
    // offsetting into the output floating-point array
    fdataptr = fdata + numbuffered * s->num_channels;
    numtoread = s->spectra_per_subint - numbuffered;
    
    // Make sure our current file number is valid 
    if (currentfile >= s->num_files)
        return 0;

    // First, attempt to read data from the current file
    numread = h5->current_sample + numtoread > s->N ? s->N - h5->current_sample : numtoread;
    vector<size_t> pos (2);
    pos[0] = h5->current_sample;
    pos[1] = 0;
    h5->bf_stokes->get2D (pos, fdata, numread, s->num_channels);
    h5->current_sample += numread;

    if (s->flip_bytes) { // byte-swap if necessary
     /* Need to add this later */
    }

    if (numread==numtoread) {  // Got all we needed
        *padding = 0;
        numbuffered = 0;
        currentblock++;
        goto return_block;
    } else {  // Didn't get all the data we needed
        numbuffered += numread;
        if (h5->current_sample >= s->N) { // End of file?
            numtopad = s->num_pad[currentfile] - numpadded;
            if (numtopad==0) {  // No padding needed.  Try reading the next file
                currentfile++;
                return get_LOFARHDF5_rawblock(fdata, s, padding);
            } else {   // Pad the data instead
                *padding = 1;
                fdataptr = fdata + numbuffered * s->num_channels;
                if (numtopad >= numtoread - numread) { // Need lots of padding
                    // Fill the rest of the buffer with padding
                    numtopad = s->spectra_per_subint - numbuffered;
                    add_padding(fdataptr, s->padvals, s->num_channels, numtopad);
                    numpadded += numtopad;
                    numbuffered = 0;
                    currentblock++;
                    // If done with padding reset padding variables
                    if (numpadded == s->num_pad[currentfile]) {
                        numpadded = 0;
                        currentfile++;
                    }
                    goto return_block;
                } else {  // Need < 1 block (or remaining block) of padding
                    add_padding(fdataptr, s->padvals, s->num_channels, numtopad);
                    numbuffered += numtopad;
                    // Done with padding, so reset padding variables
                    numpadded = 0;
                    currentfile++;
                    return get_LOFARHDF5_rawblock(fdata, s, padding);
                }
            }
        } else {
            fprintf(stderr, "Error: Problem reading record from LOFAR HDF5 data file:\n");
            fprintf(stderr, "   currentfile = %d, currentblock = %d.  Exiting.\n",
                   currentfile, currentblock);
            exit(1);
        }
    }

return_block:
    // Apply the corrections that need a full block

    // Invert the band if requested 
    if (s->apply_flipband) flip_band(fdata, s);

    // Perform Zero-DMing if requested
    if (s->remove_zerodm) remove_zerodm(fdata, s);

    return 1;
}
