import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def plotDiagnostics(dataDF, bands= ['g', 'r', 'i', 'z', 'y']): # assume cleaned up data
    ##########################################
    # plot ra, dec. in degrees.
    print 'Plotting ra, dec plot.'
    sns.jointplot(x=dataDF['ra'], y=dataDF['dec'], kind="hex", color="k")
    plt.show()
    
    ##########################################
    def plotHists(keys, df, xlabel):
        print 'Plotting %s distributions.'%xlabel
        if all(i in dataDF.keys() for i in neededKeys): # the keys needed are available
            fig, ax = plt.subplots()
            for key in neededKeys:
                plt.hist(dataDF[key], label= key, histtype= 'step', alpha= 1., lw= 3, bins= 20)
            ax.legend()
            ax.set_xlabel(xlabel)
            plt.show()
        else:
            print 'Not plotting %s since dont have the needed columns:\n%s.\n'%(xlabel, neededKeys)
        
    ##########################################
    # plot the cmodel_mag distributions
    neededKeys= ['%scmodel_mag'%b for b in bands]
    plotHists(neededKeys, dataDF, 'cmodel_mag')
    
    ##########################################
    # plot the cmodel_mag_err distributions
    neededKeys= ['%scmodel_mag_err'%b for b in bands]
    plotHists(neededKeys, dataDF, 'cmodel_mag_err')
    
    ##########################################
    # plot the absorption distributions
    neededKeys= ['a_%s'%b for b in bands]
    plotHists(neededKeys, dataDF, 'absorption')
    
    ##########################################
    # plot the <>countinputs distributions
    neededKeys= ['%scountinputs'%b for b in bands]
    plotHists(neededKeys, dataDF, 'countinputs')
  
