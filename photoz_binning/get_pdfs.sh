#!/bin/bash

cd /global/cscratch1/sd/awan/hsc_pdfs

for i in wide_aegis wide_gama09h wide_gama15h wide_hectomap wide_vvds wide_wide12h wide_xmmlss deep_cosmos deep_elaisn1 deep_xmmlss deep_deep23
do
    # create a new folder for this field
    echo starting with: $i
    #mkdir $i
    cd $i
    # get PDFs from different PZ algorithsm
    for j in nnpz #ephor ephor_ab demp frankenz
    do
        # create a new folder for this algorithm
        mkdir $j
        cd $j
        echo Getting $j/pdr1_"$j"_$i.tar.xz
        wget https://hsc-release.mtk.nao.ac.jp/archive/photoz/pdr1/pdf/$j/pdr1_"$j"_$i.tar.xz
        echo upzipping file
        tar -xf pdr1_"$j"_$i.tar.xz
        cd /global/cscratch1/sd/awan/hsc_pdfs/$i
    done
    cd /global/cscratch1/sd/awan/hsc_pdfs
done

# add group permissions
chgrp -R lsst /global/cscratch1/sd/awan/hsc_pdfs
chmod -R g-w /global/cscratch1/sd/awan/hsc_pdfs

