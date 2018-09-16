#!/bin/bash

source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh

#mkdir sl_scripts
cd sl_scripts

for i in wide_aegis wide_gama09h wide_gama15h wide_hectomap wide_vvds wide_wide12h wide_xmmlss deep_cosmos deep_elaisn1 deep_xmmlss deep_deep23
do
    for j in nnpz #ephor ephor_ab demp frankenz
    do
        cat > ${i}_${j}.sl << EOF
#!/bin/bash -l

#SBATCH --nodes=1               # Use 1 node; has multiple cores
#SBATCH --qos=regular
#SBATCH -t 00:30:00             # Set 30 min time limit
#SBATCH --constraint=haswell    # Use Haswell nodes
#SBATCH --output=/global/cscratch1/sd/awan/hsc_matched_pdfs/sbatch_output/match_pdfs.%j.out
#SBATCH --job-name=match_pdfs

srun python /global/homes/a/awan/LSST/lsstRepos/HyperSupremeStructure-HSC-LSS/photoz_binning/match_pdfs.py --fields=$i --PZalg=$j
EOF
        sbatch ${i}_${j}.sl
        echo Job submitted for ${i}_${j}
    done
done

# change permissions on the saved output. need to run these when the jobs are finished.
# chgrp -R lsst /global/cscratch1/sd/awan/hsc_matched_pdfs
# chmod -R g-w /global/cscratch1/sd/awan/hsc_matched_pdfs