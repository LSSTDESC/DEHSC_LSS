import os
import numpy as np

urls,names,tables=np.genfromtxt("urls.txt",dtype='str',unpack=True)
names=np.atleast_1d(names); urls=np.atleast_1d(urls); tables=np.atleast_1d(tables);

preurl='https://hsc-release.mtk.nao.ac.jp/datasearch/catalog_jobs/download/'

for n,u,t in zip(names,urls,tables) :
    fname_final="HSC_"+n+"_"+t+".fits"
    print(fname_final)
    if not os.path.isfile(fname_final) :
        url=preurl+u
        command1="wget "+url
        command2="mv "+u+" "+fname_final
        print(command1)
        os.system(command1)
        print(command2)
        os.system(command2)
    else :
        print("Found "+fname_final)
    os.system("mv "+fname_final+" /global/cscratch1/sd/damonge/HSC/")
