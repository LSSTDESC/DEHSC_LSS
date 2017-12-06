#!/bin/bash
export PATH=/home/lss/local/bin:$PATH
export LD_LIBRARY_PATH=/home/lss/local/lib:$LD_LIBRARY_PATH
export MPLBACKEND=PS
export PYTHONPATH=/home/lss/.local/lib/python2.7/site-packages
cd /io
python ./sample_nmt.py
