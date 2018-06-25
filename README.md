# HSC LSS analyses

We have processed the HSC data for clustering analyses following a number of steps:
1. Download the relevant data from the DR1 database. This was done using submit_job.py. This script starts all the necessary jobs in your HSC space.
   All the different fields are then downloaded using the script dwl.py with data from urls.txt.
   The full dataset is currently stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_*.fits`
2. Reduce the data. This implies:
   - Removing all the unnecessary clutter from the raw files
   - Applying sanity cuts (duplicates, etc. see description in 1705.06745).
   - Constructing maps of all relevant systematics (X-sigma depth, extinction, stars, B.O. mask).
   - Applying a magnitude cut and the star-galaxy separator.
   - Writing the reduced catalog into new files.
   All this is done with the script `process.py`. Run `python process.py -h` to see all available command-line options.
   The reduced files are stored (and available) at `/global/cscratch1/sd/damonge/HSC/HSC_processed/` in per-field directories.
3. It's worth noting that we currently use WCS to create flat-sky maps of different quantities (depth, mask, dust etc) when processing each field. The maps are generated using gnomonic projection, with the median coordinates of all sources in each field as the tangent point.
4. The bash script `run_process_all.sh` processes all the WIDE fields.
5. Once all fields have been processed, we compute maps of the galaxy distribution and their corresponding power spectra for each field. The script `study_power_spectra.py` does this for all the WIDE fields. The power spectra are contaminant-deprojected for all contaminants studied in the previous step. We expect to extend this study to a larger list of systematics.
6. Our studies currently use a magnitude limit i<24.5. This is based on a study of the 10-sigma depth maps on all the different fields, and corresponds to a conservative estimate of the magnitude limit of the sample. Note that the quality of the photo-zs degrades significantly for fainter sources (according to the HSC papers).

The scripts described above make use of some dependencies and python modules written explicitly for this work. The most relevant ones are:
- `flatmaps.py`: contains routines to construct and manipulate flat-sky maps (mimicking as much as possible the functionality implemented in HEALPix for curved skies).
- `createMaps.py`: contains routines to generate maps based on information defined on a discrete set of points.
- `dataCleanUp.py`: summarizes the steps needed to remove all useless clutter from a raw database file.
- `estDepth.py`: contains routines to create depth maps using 3 different methods. All routines are wrapped into a single one called `get_depth`.
- `flatMask.py`: describes the method used to generate the bright-object mask from the catalog data.
- `rotate.py`: contains routines to rotate a given field onto the equator.
- `NaMaster` (https://github.com/damonge/NaMaster): a python module to compute arbitrary-spin power spectra of masked fields both in flat and curved skies.

The repo currently also hosts a number of ipython notebooks that were used to carry out the first analyses on the data. These also illustrate the use that has been made of the database information.
