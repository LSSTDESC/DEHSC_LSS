#!/bin/bash

source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh

mkdir -p sl_scripts
cd sl_scripts

for i in wide_aegis wide_gama09h wide_gama15h wide_hectomap wide_vvds wide_wide12h wide_xmmlss deep_cosmos deep_elaisn1 deep_xmmlss deep_deep23
do
    for j in nnpz ephor ephor_ab demp frankenz
    do
        cat > ${i}_${j}.sl << EOF
#!/bin/bash -l

#SBATCH --nodes=1               # Use 1 node; has multiple cores
#SBATCH --qos=regular
#SBATCH -t 00:30:00             # Set 30 min limit
#SBATCH --constraint=haswell    # Use Haswell nodes
#SBATCH --output=/global/cscratch1/sd/damonge/HSC/%j_match_pdfs_${i}_${j}.out
#SBATCH --job-name=match_pdfs

srun python /global/homes/d/damonge/LSST/LSS_HSC/HyperSupremeStructure-HSC-LSS/photoz_binning/match_pdfs.py \
                    --data_main_path='/global/cscratch1/sd/damonge/HSC/HSC_processed' \
                    --pdfs_main_path='/global/cscratch1/sd/awan/hsc_pdfs/' \
                    --outDir='/global/cscratch1/sd/damonge/HSC/HSC_processed/${i^^}' \
                    --fields=$i --PZalg=$j
EOF
        sbatch ${i}_${j}.sl
        echo Job submitted for ${i}_${j}
    done
done

# change permissions on the saved output. need to run these when the jobs are finished.
# chgrp -R lsst /global/cscratch1/sd/awan/hsc_matched_pdfs
# chmod -R g-w /global/cscratch1/sd/awan/hsc_matched_pdfs
