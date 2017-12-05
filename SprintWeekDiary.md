# Tasks carried out during the 2017 sprint week

## Day 1
We downloaded the catalogs for the GAMA15H field. 

Camille, Humna and Javi familiarized themselves with the data, making the first awesome plots. Results SampleDiagnosticsForHSCDR1.ipynb (and Humna's forked version of the same file)

## Day 2
Both the forced and random catalogs can be found in /global/cscratch1/sd/damonge/HSC/

We have decided to work with flat-sky maps, in order to speed up the progress. David will work on building the framework to deal with flat-sky maps. For the moment a skeleton structure can be found in flatmaps.py (see sample_fm.py for an example of how to use it).

In the meantime, Camille, Humna and Javi will work on generating a set of interesting maps from the data and randoms.
Good candidates for these maps would be
* Depth (limiting magnitude). It'd be good to check several versions of this map:
  * Using Javi's method
  * Binning 5*flux_error for all galaxies in each pixel. In this case, we'd like to look at the plot of this quantity and of its variance.
  * Using the random sky_std as described in https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points-for-dr1/. It'd be great to have versions of this map for all filters
* PSF. This can be done from the randoms (there's an ipsf_size column)
* Xcountinputs (where X is g,r,i,z,y) -> number of images contributing to each object
* Dust extinction (a_X in the forced table)
* Star map (from iclassification_extendedness=0)
* Bright object masks: check out the flags iflags_pixel_bright_object_center and iflags_pixel_bright_object_any, either in the random or forced catalogs (actually, it would be nice that both catalogs produce the same mask).
* Airmass: we don't have a clear idea of how to get that map yet, so any ideas will be welcome!

 