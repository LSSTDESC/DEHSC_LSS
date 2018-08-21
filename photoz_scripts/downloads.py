import urllib
import gzip
import os


def get_COSMOS_photoz_cat(url = 'ftp://ftp.iap.fr/pub/from_users/hjmcc/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits.gz'):

    # Downloads and extracts the COSMOS photoz catalog

    print 'Downloading COSMOS2015_Laigle+_v1.1.fits.gz...'
    urllib.urlretrieve(url, 'COSMOS2015_Laigle+_v1.1.fits.gz')

    print 'Decompressing COSMOS2015_Laigle+_v1.1.fits.gz...'
    with gzip.open('./COSMOS2015_Laigle+_v1.1.fits.gz', 'rb') as readfile:
        with open('./COSMOS2015_Laigle+_v1.1.fits', 'wb') as writefile:
            gzdata = readfile.read()
            writefile.write(gzdata)

    os.remove('./COSMOS2015_Laigle+_v1.1.fits.gz')


def get_HSC_pdfs(url = 'https://hsc-release.mtk.nao.ac.jp/archive/photoz/pdr1/pdf/ephor/pdr1_ephor_deep_cosmos.tar.xz'):

    # Downloads and extracts the HSC pdf data

    print 'Downloading ' + url.split('/')[-2] + ' PDFs...'
    urllib.urlretrieve(url, 'pdr1_ephor_deep_cosmos.tar.xz')

    print 'Decompressing pdr1_ephor_deep_cosmos.tar.xz...'
    os.system('tar -xJf pdr1_ephor_deep_cosmos.tar.xz')

    os.remove('./pdr1_ephor_deep_cosmos.tar.xz')