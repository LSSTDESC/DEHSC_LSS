from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from glob import glob
from astropy.io import fits
from tqdm import tqdm



class reader:

    def __init__(self, inputfp = './pdr1_ephor_deep_cosmos/'):

        file_list = sorted(glob(inputfp + '*'))

        self.object_id = np.array([])
        self.pdf = []
        self.bins = np.array([])

        idlist = []
        pdflist = []
        binlist = []

        for filename in tqdm(file_list):

            hdulist = fits.open(filename)

            idlist = idlist + list(hdulist[1].data['ID'])
            pdflist.append(hdulist[1].data['PDF'])
            binlist.append(hdulist[2].data['BINS'])

        if all([all(binlist[0] == rest) for rest in binlist]):
            self.bins = binlist[0]

        self.pdf = np.concatenate(tuple(pdflist))

        self.object_id = np.array(idlist)