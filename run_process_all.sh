#!/bin/bash

predir_out=/global/cscratch1/sd/damonge/HSC

#module swap PrgEnv-gnu PrgEnv-intel
#module unload python/3.6-anaconda-4.4
#module load python/2.7-anaconda-4.4

#First clean up the metadata
for table in WIDE DEEP UDEEP
do
    python process_metadata.py --input-file ${predir_out}/HSC_${table}_frames.fits --output-file ${predir_out}/HSC_processed/HSC_${table}_frames_proc.fits
done

#Now clean up the WIDE fields
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    mkdir -p ${dirname}
    #python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --gen-plots --min-snr 10.0 --depth-cut 24.5
    python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --min-snr 10.0 --depth-cut 24.5 --flat-project CAR
done

#Same for deep fields.
for field in DEEP_COSMOS DEEP_DEEP23 DEEP_ELAISN1 DEEP_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    mkdir -p ${dirname}
    #python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --gen-plots --min-snr 10.0 --depth-cut 24.5
    python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --min-snr 10.0 --depth-cut 24.5 --flat-project CAR
done
exit

#Compute per-frame systematic maps
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    echo $field
    dirname=${predir_out}/HSC_processed/${field}
    python map_obscond.py --input-frames ${predir_out}/HSC_processed/HSC_WIDE_frames_proc.fits --map-sample ${dirname}/${field}_MaskedFraction.fits --output-prefix ${dirname}/${field}
done

#Compute galaxy count maps for each field
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    echo $field
    dirname=${predir_out}/HSC_processed/${field}
    python cat_sampler.py --input-prefix ${dirname}/${field} --output-file ${dirname}/${field}_Ngal_bins_eab_best_single.fits --pz-type ephor_ab --pz-mark best --pz-bins bins_z_single.txt --map-sample ${dirname}/${field}_MaskedFraction.fits --analysis-band i --depth-cut 24.5
done


#Compute cross-power spectra for each field
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    python power_specter.py --output-file ${dirname}/${field}_spectra_eab_best_single_cont_dpt_dst_str_ams_fwh_ssk.sacc --input-prefix ${dirname}/${field} --input-maps ${dirname}/${field}_Ngal_bins_eab_best_single.fits --ell-bins ell_bins.txt --mcm-output ${dirname}/${field}_mcm.dat --hsc-field HSC_${field} --cont-depth --cont-dust --cont-stars --cont-oc airmass,seeing,sigma_sky --covariance-option theory --covariance-theory-prediction cl_th_single.txt --covariance-coupling-file ${dirname}/${field}_cov_mcm.dat
    python power_specter.py --output-file ${dirname}/${field}_spectra_eab_best_single_nocont.sacc --input-prefix ${dirname}/${field} --input-maps ${dirname}/${field}_Ngal_bins_eab_best_single.fits --ell-bins ell_bins.txt --mcm-output ${dirname}/${field}_mcm.dat --hsc-field HSC_${field} --covariance-option theory --covariance-theory-prediction cl_th_single.txt --covariance-coupling-file ${dirname}/${field}_cov_mcm.dat
done

#So far we've only looked at the WIDE fields
exit

for field in COSMOS_WIDE_BEST COSMOS_WIDE_MEDIAN COSMOS_WIDE_WORST 
do
    dirname=${predir_out}/HSC_processed/${field}
    mkdir -p ${dirname}
    python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --gen-plots --min-snr 10.0 --depth-cut 24.5
done

for field in DEEP_COSMOS DEEP_DEEP32 DEEP_ELAISN1 DEEP_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    mkdir -p ${dirname}
    python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --gen-plots --min-snr 10.0 --depth-cut 24.5
done

for field in UDEEP_COSMOS UDEEP_SXDS
do
    dirname=${predir_out}/HSC_processed/${field}
    mkdir -p ${dirname}
    python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --gen-plots --min-snr 10.0 --depth-cut 24.5
done
