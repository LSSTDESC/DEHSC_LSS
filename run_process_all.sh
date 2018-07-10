#!/bin/bash

predir_out=/global/cscratch1/sd/damonge/HSC

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
    python process.py --input-field ${field} --resolution 0.01 --field-padding 0.1 --output-prefix ${dirname}/${field} --save-systematics --save-masks --save-depth-maps --min-snr 10.0 --depth-cut 24.5
done

#Compute per-frame systematic maps
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    python map_obscond.py --input-frames ${predir_out}/HSC_processed/HSC_WIDE_frames_proc.fits --map-sample ${dirname}/${field}_MaskedFraction.fits --output-prefix ${dirname}/${field}
done

#Compute galaxy count maps for each field
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    python cat_sampler.py --input-prefix ${dirname}/${field} --output-file ${dirname}/${field}_Ngal_bins_eab_best_two.fits --pz-type ephor_ab --pz-mark best --pz-bins bins_z_two.txt --map-sample ${dirname}/${field}_MaskedFraction.fits --analysis-band i --depth-cut 24.5
done

#Compute cross-power spectra for each field
for field in WIDE_AEGIS WIDE_GAMA09H WIDE_GAMA15H WIDE_HECTOMAP WIDE_VVDS WIDE_WIDE12H WIDE_XMMLSS
do
    dirname=${predir_out}/HSC_processed/${field}
    echo python power_specter.py --output-file ${dirname}/${field}_spectra_eab_best_two.sacc --input-prefix ${dirname}/${field} --input-maps ${dirname}/${field}_Ngal_bins_eab_best_two.fits --ell-bins ell_bins.txt --mcm-output ${dirname}/${field} --hsc-field HSC_${field} --cont-depth --cont-dust --cont-stars --cont-oc airmass,seeing,sigma_sky
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
