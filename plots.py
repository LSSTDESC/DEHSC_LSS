from __future__ import print_function
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
import os

def plotDiagnostics(dataDF, bands= ['g', 'r', 'i', 'z', 'y']): # assume cleaned up data
    ##########################################
    # plot ra, dec. in degrees.
    print('Plotting ra, dec plot.')
    sns.jointplot(x=dataDF['ra'], y=dataDF['dec'], kind="hex", color="k")
    plt.show()
    
    ##########################################
    def plotHists(keys, df, xlabel):
        print('Plotting %s distributions.'%xlabel)
        if all(i in dataDF.keys() for i in neededKeys): # the keys needed are available
            fig, ax = plt.subplots()
            for key in neededKeys:
                plt.hist(dataDF[key], label= key, histtype= 'step', alpha= 1., lw= 3, bins= 20)
            ax.legend()
            ax.set_xlabel(xlabel)
            plt.show()
        else:
            print('Not plotting %s since dont have the needed columns:\n%s.\n'%(xlabel, neededKeys))
        
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
  
def plot_wtheta(theta, wtheta, wtheta_sig, title= None):
    # theta assumed to be in degrees
    if title is None: title= ''

    plt.clf()
    fig, ax= plt.subplots(1,1)
    ax.plot(theta, np.zeros(len(theta))+1, color= 'k', lw=2, linestyle= ':')
    ax.plot(theta, np.zeros(len(theta)), color= 'k', lw=2, linestyle= ':')
    ax.plot(theta, np.zeros(len(theta))-1, color= 'k', lw=2, linestyle= ':')

    color= 'b'
    ax.plot(theta, wtheta, lw= 2, color= color)
    ax.scatter(theta, wtheta, color= color)
    ax.errorbar(theta, wtheta, yerr= wtheta_sig, color=color)

    fontsize= 12
    ax.set_title(title, fontsize= fontsize)
    ax.set_xlabel(r'$\theta$ (Degrees)', fontsize= fontsize+2)
    ax.set_ylabel(r'$w(\theta)$', fontsize= fontsize+2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    plt.show()
    
def plot_wthetas(theta_dict, wtheta_dict, wtheta_sig_dict, title= None,
                 savePlot= False, outDir= None, filename= None):
    # theta assumed to be in degrees
    if title is None: title= ''

    colors= ['b', 'r', 'g', 'k', 'm', 'c']
    
    plt.clf()
    fig, ax= plt.subplots(1,1)
    for i, tag in enumerate(theta_dict):
        ax.plot(theta_dict[tag], np.zeros(len(theta_dict[tag]))+1, color= 'k', lw=2, linestyle= ':')
        ax.plot(theta_dict[tag], np.zeros(len(theta_dict[tag])), color= 'k', lw=2, linestyle= ':')
        ax.plot(theta_dict[tag], np.zeros(len(theta_dict[tag]))-1, color= 'k', lw=2, linestyle= ':')

        if (i>len(colors)-1):
            style= '.-'
            i= i%len(colors)
        else:
            style= 'o-'
        ax.plot(theta_dict[tag], wtheta_dict[tag], style, lw= 2, color= colors[i], label= tag)
        ax.errorbar(theta_dict[tag], wtheta_dict[tag], yerr= wtheta_sig_dict[tag], color=colors[i])

    fontsize= 12
    ax.legend(fontsize= fontsize, bbox_to_anchor= (1.4, 0.5))
    ax.set_title(title, fontsize= fontsize)
    ax.set_xlabel(r'$\theta$ (Degrees)', fontsize= fontsize+2)
    ax.set_ylabel(r'$w(\theta)$', fontsize= fontsize+2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    
    if savePlot:
        if outDir is None:
            raise ValueError('Need outDir if want to save plot.')
        if filename is None:
            raise ValueError('Need filename if want to save plot.')
            
        current= os.getcwd()
        os.chdir(outDir)
        plt.savefig(filename,bbox_inches='tight')
        plt.close('all')
        os.chdir(current)
        print('Saved ', filename)
    else:
        plt.show()
