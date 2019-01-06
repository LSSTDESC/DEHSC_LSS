import urllib
import gzip
import os

url = 'ftp://ftp.iap.fr/pub/from_users/hjmcc/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits.gz'

print 'Downloading COSMOS2015_Laigle+_v1.1.fits.gz...'
urllib.urlretrieve(url, 'COSMOS2015_Laigle+_v1.1.fits.gz')

print 'Decompressing COSMOS2015_Laigle+_v1.1.fits.gz...'
with gzip.open('./COSMOS2015_Laigle+_v1.1.fits.gz', 'rb') as readfile:
    with open('./COSMOS2015_Laigle+_v1.1.fits', 'wb') as writefile:
        gzdata = readfile.read()
        writefile.write(gzdata)

os.remove('./COSMOS2015_Laigle+_v1.1.fits.gz')
