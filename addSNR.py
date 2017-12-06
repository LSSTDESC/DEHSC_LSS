import matplotlib.pyplot as plt
def addSNR(dataDF, bands= ['g', 'r', 'i', 'z', 'y']):
    # calculate the SNR as cmodel_flux/flux_err for each band. add as a column to the dataframe.
    SNCols= []
    for band in bands:
        key= '%scmodel'%band
        # check is cmodel is present
        if '%s_flux'%key not in dataDF.keys():
            raise ValueError('%s not in the dataframe.'%key)
        # check is cmodel is present
        if '%s_flux_err'%key not in dataDF.keys():
            raise ValueError('%s not in the dataframe.'%key)
            
        # check is SNR is present already
        if '%s-SNR'%key in dataDF.keys():
            print '%s in the dataframe the datafram already. Overwriting.'%key
            
        SNkey= '%s-SNR'%key
        SNCols.append(SNkey)
        # if needed column is there, calculate SNR
        dataDF[SNkey]= dataDF['%s_flux'%key]/dataDF['%s_flux_err'%key]
    
    print 'Added the columns. Now plotting the distribution.'
    # plot the distbriution of the S/N
    fig, ax = plt.subplots()
    for key in SNCols:
        plt.hist(dataDF[key], label= key, histtype= 'step', alpha= 1., lw= 3,
                 bins= [0., 1, 2., 3., 4., 5., 10, 50, 500, 1000])
    ax.legend()
    ax.set_xlabel('SNR')
    ax.set_ylabel('Object count')
    ax.set_xscale('log')
    
    return dataDF, SNCols


