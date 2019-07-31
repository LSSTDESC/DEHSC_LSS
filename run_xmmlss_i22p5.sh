#!/bin/bash
  
#python3 -m hsc_lss ReduceCat   --raw_data=./hsc_lss_params/input_list_xmmlss.txt   --clean_catalog=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/clean_catalog.fits   --dust_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/dust_map.fits   --star_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/star_map.fits   --bo_mask=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/bo_mask.fits   --masked_fraction=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/masked_fraction.fits   --depth_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/depth_map.fits   --config=./hsc_lss_params/params_i22p5/config.yml 

#python3 -m hsc_lss SystMapper   --frames_data=/global/cscratch1/sd/damonge/HSC_ceci/PDR1_WIDE_frames.fits   --masked_fraction=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/masked_fraction.fits   --ccdtemp_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ccdtemp_maps.fits   --airmass_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/airmass_maps.fits   --exptime_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/exptime_maps.fits   --skylevel_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/skylevel_maps.fits   --sigma_sky_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/sigma_sky_maps.fits   --seeing_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/seeing_maps.fits   --ellipt_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ellipt_maps.fits   --nvisit_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/nvisit_maps.fits   --config=./hsc_lss_params/params_i22p5/config.yml 

#python3 -m hsc_lss PDFMatch   --clean_catalog=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/clean_catalog.fits   --pdf_dir=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS   --pdf_matched=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/pdf_matched.txt   --config=./hsc_lss_params/params_i22p5/config.yml 

#python3 -m hsc_lss COSMOSWeight   --cosmos_data=/global/cscratch1/sd/damonge/HSC_ceci/COSMOS2015_Laigle+_v1.1.fits   --cosmos_hsc=/global/cscratch1/sd/damonge/HSC_ceci/PDR1_DEEP_COSMOS_forced.fits   --cosmos_weights=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/cosmos_weights.fits   --config=./hsc_lss_params/params_i22p5/config.yml 

#python3 -m hsc_lss CatMapper   --clean_catalog=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/clean_catalog.fits   --masked_fraction=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/masked_fraction.fits   --cosmos_weights=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/cosmos_weights.fits   --pdf_matched=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/pdf_matched.txt   --ngal_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ngal_maps.fits   --config=./hsc_lss_params/params_i22p5/config.yml 

#python3 -m hsc_lss MapDiagnoser   --masked_fraction=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/masked_fraction.fits   --ngal_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ngal_maps.fits   --dust_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/dust_map.fits   --star_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/star_map.fits   --depth_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/depth_map.fits   --ccdtemp_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ccdtemp_maps.fits   --airmass_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/airmass_maps.fits   --exptime_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/exptime_maps.fits   --skylevel_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/skylevel_maps.fits   --sigma_sky_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/sigma_sky_maps.fits   --seeing_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/seeing_maps.fits   --ellipt_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ellipt_maps.fits   --nvisit_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/nvisit_maps.fits   --systmap_plots=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/systmap_plots   --config=./hsc_lss_params/params_i22p5/config.yml 

#python3 -m hsc_lss PowerSpecter   --masked_fraction=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/masked_fraction.fits   --ngal_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ngal_maps.fits   --dust_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/dust_map.fits   --star_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/star_map.fits   --depth_map=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/depth_map.fits   --ccdtemp_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ccdtemp_maps.fits   --airmass_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/airmass_maps.fits   --exptime_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/exptime_maps.fits   --skylevel_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/skylevel_maps.fits   --sigma_sky_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/sigma_sky_maps.fits   --seeing_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/seeing_maps.fits   --ellipt_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/ellipt_maps.fits   --nvisit_maps=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/nvisit_maps.fits   --cosmos_weights=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/cosmos_weights.fits   --syst_masking_file=./hsc_lss_params/systematic_cuts/WIDE_XMMLSS_syst_cuts.txt   --dummy=/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i22p5_out/dummy   --config=./hsc_lss_params/params_i22p5/config.yml 

