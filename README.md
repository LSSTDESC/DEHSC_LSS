# HSC LSS analyses

We have processed the HSC data for clustering analyses following a number of steps:
1. Download the relevant data from the DR1 database. This was done using the scripts in the directory sql_utils. The script submit_job.py starts all the necessary jobs in your HSC space.
   All the different fields are then downloaded using the script dwl.py with data from urls.txt.
   The full dataset is currently stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_*.fits`. This includes both the forced photometry catalog and the metadata.
2. Reduce the metadata. This implies:
   - Removing all the unnecessary clutter from the raw files
   - Filter out all unnecessary columns (see line 50 of `process_metadata.py` for the columns we actually keep).
   - Writing the reduced tables into new files.
   All this is done with the script `process_metadata.py`. Run `python process_metadata.py -h` to see all available command-line options. The reduced files are stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_processed/HSC_*_frames_proc.fits`.
3. Reduce the catalog data. This implies:
   - Removing all the unnecessary clutter from the raw files
   - Applying sanity cuts (duplicates, etc. see description in 1705.06745).
   - Constructing maps of all relevant catalog-based systematics (X-sigma depth, extinction, stars, B.O. mask).
   - Applying a magnitude cut and the star-galaxy separator.
   - Writing the reduced catalog into new files.
   All this is done with the script `process.py`. Run `python process.py -h` to see all available command-line options.
   The reduced files are stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_processed/` in per-field sub-directories.
4. Use the metadata to generate maps of the per-frame systematics (i.e. observing conditions) in each field. This implies:
   - Selecting the frames that fall within the field region.
   - Computing the map pixels each frame intersects and their corresponding areas (this is the most time-consuming part).
   - For each map pixel, build a histogram of the values of each relevant systematic in each exposure touching that pixel.
   - Compress those histograms into maps of summary statistics (currently only computing the mean of each quantity).
   - Maps are generated per-band.
   - The currently mapped quantities are given in line 13 of `map_obscond.py`.
   All this is done with the script `map_obscond.py`. Run `python map_obscond.py -h` to see all available command-line options. The maps are stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_processed/` in the sub-directories created above for each field.
5. Create galaxy count maps in a number of redshift bins for each field. This implies:
   - Reading in the processed catalogs.
   - Binning them in terms of a given photo-z marker (e.g. ML, mean etc.) for a particular photo-z code.
   - Generating a map of the number of objects in a given bin found per pixel.
   - Generating an estimate of the redshift distribution for objects in a given bin. We currently do this as a histogram of the MC redshift value stored for each object.
   - Save all maps and redshift distributions to file. These are currently collected into a single FITS file that alternates image HDUs (containing the maps) and table HDUs (containing the binned N(z)).
   All this is done with the script `cat_sampler.py`. Run `python cat_sampler.py -h` to see all possible command-line options. The results are stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_processed/` in the sub-directories created above for each field.
6. It's worth noting that we currently use WCS to create flat-sky maps of different quantities (depth, mask, dust etc) when processing each field. The maps are generated using gnomonic projection, with the median coordinates of all sources in each field as the tangent point.
7. The bash script `run_process_all.sh` runs 2, 3, 4 and 5 above for all the WIDE fields.
8. Once all fields have been processed, we compute maps of the galaxy distribution and their corresponding power spectra for each field. The script `study_power_spectra.py` does this for all the WIDE fields. The power spectra are contaminant-deprojected for all contaminants studied in the previous step. We expect to extend this study to a larger list of systematics.
9. Our studies currently use a magnitude limit i<24.5. This is based on a study of the 10-sigma depth maps on all the different fields, and corresponds to a conservative estimate of the magnitude limit of the sample. Note that the quality of the photo-zs degrades significantly for fainter sources (according to the HSC papers).

The scripts described above make use of some dependencies and python modules written explicitly for this work. The most relevant ones are:
- `flatmaps.py`: contains routines to construct and manipulate flat-sky maps (mimicking as much as possible the functionality implemented in HEALPix for curved skies).
- `createMaps.py`: contains routines to generate maps based on information defined on a discrete set of points.
- `dataCleanUp.py`: summarizes the steps needed to remove all useless clutter from a raw database file.
- `estDepth.py`: contains routines to create depth maps using 3 different methods. All routines are wrapped into a single one called `get_depth`.
- `flatMask.py`: describes the method used to generate the bright-object mask from the catalog data.
- `rotate.py`: contains routines to rotate a given field onto the equator.
- `NaMaster` (https://github.com/damonge/NaMaster): a python module to compute arbitrary-spin power spectra of masked fields both in flat and curved skies.

The repo currently also hosts a number of ipython notebooks that were used to carry out the first analyses on the data. These also illustrate the use that has been made of the database information.
