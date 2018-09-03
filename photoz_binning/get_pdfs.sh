#!/bin/bash

cd /global/cscratch1/sd/awan/hsc_pdfs

for i in aegis gama09h gama15h hectomap vvds wide12h xmmlss
do
        echo starting with: $i
        mkdir $i
        cd $i
        wget https://hsc-release.mtk.nao.ac.jp/archive/photoz/pdr1/pdf/ephor_ab/pdr1_ephor_ab_wide_$i.tar.xz
        echo upzipping file
        tar -xf pdr1_ephor_ab_wide_$i.tar.xz
        
        #wget https://hsc-release.mtk.nao.ac.jp/archive/photoz/pdr1/pdf/nnpz/pdr1_nnpz_wide_$i.tar.xz
        #tar -xf pdr1_nnpz_wide_$i.tar.xz
        
        #wget https://hsc-release.mtk.nao.ac.jp/archive/photoz/pdr1/pdf/frankenz/pdr1_frankenz_wide_$i.tar.xz
        #tar -xf pdr1_frankenz_wide_$i.tar.xz
        
        cd /global/cscratch1/sd/awan/hsc_pdfs
done



