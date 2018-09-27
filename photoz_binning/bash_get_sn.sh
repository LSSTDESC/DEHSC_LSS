#!/bin/bash

source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh

mkdir -p sl_scripts
cd sl_scripts

#########################################################################################################
# run the analysis on the three wide fields using dn/dz from z_mc
for i in wide_aegis wide_vvds wide_xmmlss
do
    for j in best mode
    do
        cat > ${i}_z${j}_nz-mc.sl << EOF
#!/bin/bash -l

#SBATCH --nodes=1               # Use 1 node; has multiple cores
#SBATCH --qos=regular
#SBATCH -t 00:30:00             # Set 30 min time limit
#SBATCH --constraint=haswell    # Use Haswell nodes
#SBATCH --output=/global/cscratch1/sd/awan/lsst_output/hsc_output/sbatch_output/%j_get_photoz_sn_${i}_z${j}_nz-mc.out
#SBATCH --job-name=pz_sn

srun python /global/homes/a/awan/LSST/lsstRepos/HyperSupremeStructure-HSC-LSS/photoz_binning/get_sn_photoz_bins.py \
                                    --fields=$i --nbin=10 --z_type=$j \
                                    --dont_show_plots --save_plots --nz_mc \
                                    --outDir='/global/cscratch1/sd/awan/lsst_output/hsc_output'
EOF
        sbatch ${i}_z${j}_nz-mc.sl
        echo Job submitted for ${i}_z${j}_nz-mc
    done
done

#########################################################################################################
# for comparison, run the analysis on wide_aegis using dn/dz from pdf stacking
for i in wide_aegis
do
    for j in best mode
    do
        cat > ${i}_z${j}_nz-pdfs.sl << EOF
#!/bin/bash -l

#SBATCH --nodes=1               # Use 1 node; has multiple cores
#SBATCH --qos=regular
#SBATCH -t 00:30:00             # Set 30 min time limit
#SBATCH --constraint=haswell    # Use Haswell nodes
#SBATCH --output=/global/cscratch1/sd/awan/lsst_output/hsc_output/sbatch_output/%j_get_photoz_sn_${i}_z${j}_nz-mc.out
#SBATCH --job-name=pz_sn

srun python /global/homes/a/awan/LSST/lsstRepos/HyperSupremeStructure-HSC-LSS/photoz_binning/get_sn_photoz_bins.py \
                                    --fields=$i --nbin=10 --z_type=$j \
                                    --dont_show_plots --save_plots \
                                    --outDir='/global/cscratch1/sd/awan/lsst_output/hsc_output'
EOF
        sbatch ${i}_z${j}_nz-pdfs.sl
        echo Job submitted for ${i}_z${j}_nz-pdfs
    done
done

# change permissions on the saved output. need to run these when the jobs are finished.
# chgrp -R lsst /global/cscratch1/sd/awan/hsc_matched_pdfs
# chmod -R g-w /global/cscratch1/sd/awan/hsc_matched_pdfs