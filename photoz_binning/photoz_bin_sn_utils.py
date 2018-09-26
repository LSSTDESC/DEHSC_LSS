import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt

__all__ = ['get_dn_dz', 'get_bins_ends', 'get_bin_edges',
          'get_surface_number_density', 'calc_sn']

# --------------------------------------------------------------------------------------------------------
def get_dn_dz(z_edges, hsc_z_phot, hsc_ids, matched_pdf_ids, matched_pdfs, nz_mc=False, hsc_z_mc=None, z_bins=None):
    # to get dn/dz, we first find all the galaxies whose z_phot are in a given bin
    # then dn/dz is the sum of the pdfs of these galaxies or a histogram of their z_mc
    print('Running get_dn_dz ... ')
    if nz_mc:
        if hsc_z_mc is None:
            raise ValueError('Need z_mc to estimate dn/dz.')
        if z_bins is None:
            raise ValueError('Need the bins array to estimates dn/dz.')
        # find the bin centers
        diff = np.unique([round(z_bins[i+1]-z_bins[i],2) for i in range(len(z_bins)-1)])
        if len(diff)>1:
            print('Finding multiple $\Delta$z for binning: %s. Using %s '%(diff, diff[0]))
        hist_bins = list(z_bins-diff[0]) + [max(z_bins)+diff[0]]

    # calculate dn/dz now
    dndz = {}
    num_gals = {}
    for i in range(len(z_edges)-1):
        zmin, zmax = z_edges[i], z_edges[i+1]
        bin_key = '%s\leq z<%s'%(zmin, zmax)
        dndz[bin_key] = []
        ind = np.where((hsc_z_phot >= zmin) & (hsc_z_phot < zmax))[0]
        num_gals[bin_key] = len(ind)
        # decide on the way to get dn/dz
        if nz_mc:
            print('Creating the z_mc histogram of %s objects for %s'%(num_gals[bin_key], bin_key))
            dndz[bin_key], _, _ = plt.hist(hsc_z_mc[ind], bins=hist_bins)
        else:
            print('Summing pdfs of %s objects for %s'%(num_gals[bin_key], bin_key))
            for j, ID in enumerate(hsc_ids[ind]):
                row = np.where(matched_pdf_ids==ID)[0][0]
                if j==0:
                    dndz[bin_key] = matched_pdfs[row]
                else:
                    dndz[bin_key] += matched_pdfs[row]
    return dndz, num_gals

# --------------------------------------------------------------------------------------------------------
def get_bins_ends(z_bin_array):
    # get the labels and bin edges based on the input z-bin array
    bin_keys, bin_ends = [], {}
    for i in range(len(z_bin_array)-1):
        zmin, zmax = z_bin_array[i], z_bin_array[i+1]
        bin_key = '%s\leq z<%s'%(zmin, zmax)
        bin_keys.append(bin_key)
        bin_ends[bin_key] = [zmin, zmax]
    return np.array(bin_keys), bin_ends

# --------------------------------------------------------------------------------------------------------
def get_bin_edges(nbin, hsc_z_phot, z_bins):
    # decide on bins edges: we look for bin edges that have roughly the same number of galaxies
    zmin, zmax = 0.15, 1.5   # these are fixed since this is the range where z_phot is reliable
    z_bins_finer = np.arange(min(z_bins), max(z_bins), 0.001)
    obj_thres = 50    # the bins would have galaxies within this threshold of each other
    print('obj_thres: %s'%obj_thres)
    if nbin==0:
        print('nbin must be greater than one.')
        return
    else:
        n_obj = len(np.where((hsc_z_phot >= zmin) & (hsc_z_phot < zmax))[0])
        if nbin==1:
            print('%s objects in %s<=z<%s'%(n_obj, zmin, zmax))
            bin_ends = np.array([zmin, zmax])
        else:
            # want the same number of objects in the bins
            wanted_n_obj_in_bin = n_obj/nbin
            print('Looking for bins with about %s galaxies each ... '%wanted_n_obj_in_bin)

            bin_ends = np.zeros(nbin+1)
            bin_ends[0] = zmin
            bin_ends[-1] = zmax

            for i in range(nbin-1):
                print('\nLooking for %sth bin egde'%(i+1))
                n_obj = 0
                j = 0
                good_edges = [f for f in z_bins_finer if f>bin_ends[i]]
                # loop over the bins until have the target number of galaxies
                while ((n_obj-wanted_n_obj_in_bin) < obj_thres) and (j<len(good_edges)):
                    z_edge = good_edges[j]
                    n_obj = len(np.where((hsc_z_phot >= bin_ends[i]) & (hsc_z_phot < z_edge))[0])
                    j += 1
                # some checks
                if (z_edge==bin_ends[-1]):
                    raise ValueError('Need to change the threshold on the number of galaxies: %s currently'%(obj_thres))
                if n_obj==0:
                    raise ValueError('Something is wrong. No objects found in this bin.')
                # things look ok. store.
                bin_ends[i+1] = float('%.2f'%z_edge)

            print('\nFinal:')
            for i in range(nbin):
                n_obj = len(np.where((hsc_z_phot >= bin_ends[i]) & (hsc_z_phot < bin_ends[i+1]))[0])
                print('%s objects in %s<=z<%s'%(n_obj, bin_ends[i], bin_ends[i+1]))

        return bin_ends
# --------------------------------------------------------------------------------------------------------
def get_surface_number_density(n_objs, area_in_sr):  # in 1/Sr
    return n_objs/area_in_sr

# --------------------------------------------------------------------------------------------------------
def calc_sn(z_phot, z_bins, hsc_z_phot, hsc_ids, matched_pdf_ids, matched_pdfs, n_z, ell, area_in_sr,
            plot_cls=True, hsc_z_mc=None, nz_mc=False,
            save_plots=False, dont_show_plots=False, filetag=None, outDir=None):
    print('-----------------------------------------\nRunning calc_sn for %s bins ... '%(len(z_phot)-1))
    if nz_mc:
        if hsc_z_mc is None:
            raise ValueError('Need either the N(z) or z_mc to estimate dn/dz.')

    # set up the bins
    z_edges_all = [min(z_bins)] + list(z_phot) + [max(z_bins)]   # cover the entire range, even outside the target bins

    target_bins, target_bin_ends = get_bins_ends(z_bin_array=z_phot)
    all_bins, all_bin_ends = get_bins_ends(z_bin_array=z_edges_all)

    print('Target bins: %s'%(target_bins))
    print('All bins: %s\n'%(all_bins))

    # get dn/dz
    dn_dz, num_gals = get_dn_dz(z_edges=z_edges_all, hsc_z_phot=hsc_z_phot,
                                hsc_ids=hsc_ids, matched_pdf_ids=matched_pdf_ids, matched_pdfs=matched_pdfs,
                                nz_mc=nz_mc, hsc_z_mc=hsc_z_mc, z_bins=z_bins)
    # plot the dn/dz for the different bins
    fontsize = 20
    plt.clf()
    for key in dn_dz:
        if key in target_bin_ends:
            plt.plot(z_bins, dn_dz[key], lw=4, label=r'** dn_dz: $%s$: %s galaxies'%(key, num_gals[key]))
        else:
            plt.plot(z_bins, dn_dz[key], lw=4, label=r'   dn_dz: $%s$: %s galaxies'%(key, num_gals[key]))
    # plot n(z) for reference
    plt.plot(z_bins, n_z, 'k:', label='   n_z', alpha=0.7)
    # plot bin edges
    ymin, ymax = plt.gca().get_ylim()
    for key in target_bin_ends:
        for z_edge in target_bin_ends[key]:
            plt.plot([z_edge, z_edge], [0, ymax], 'k-.', alpha=0.7)
    plt.title(filetag, fontsize=fontsize)
    plt.gca().tick_params(axis='both', labelsize=fontsize-2)
    plt.xlabel('z', fontsize=fontsize)
    plt.gcf().set_size_inches(10, 6)
    plt.legend(bbox_to_anchor=(1,1), fontsize=fontsize-5)
    if save_plots:
        filename = '%s_dndz_%sbins.png'%(filetag, len(z_phot)-1)
        plt.savefig('%s/%s'%(outDir, filename), format='png', bbox_inches='tight')
        print('\n## Saved plot: %s\n'%filename)
    if dont_show_plots:
        plt.close('all')
    else:
        plt.show()

    # set up ccl
    cosmo_fid= ccl.Cosmology(Omega_c=0.25, Omega_b=0.05, h=0.7, sigma8=0.8, n_s=0.96)
    bias = 1 + 0.84*z_bins

    # calculate the signal matrix
    print('Calculating the signal matrix ... ')
    S = np.zeros(shape=(len(ell), len(all_bins), len(all_bins)))
    for i, bin1_key in enumerate(all_bins):
        nc_bin1 = ccl.ClTracerNumberCounts(cosmo_fid, has_rsd=False, has_magnification=False,
                                            n=dn_dz[bin1_key], bias=bias, z=z_bins)
        for j, bin2_key in enumerate(all_bins):
            nc_bin2= ccl.ClTracerNumberCounts(cosmo_fid, has_rsd=False, has_magnification=False,
                                                  n=dn_dz[bin2_key], bias=bias, z=z_bins)
            S[:, i, j] = ccl.angular_cl(cosmo_fid, nc_bin1, nc_bin2, ell)

    # calculate the noise matrix
    print('Calculating the noise matrix ... ')
    N = np.zeros(shape=(len(ell), len(all_bins), len(all_bins)))
    for i, bin_key in enumerate(all_bins):
        zmin, zmax = all_bin_ends[bin_key]
        n_obj = len(np.where((hsc_z_phot > zmin) & (hsc_z_phot < zmax))[0])
        print('%s objects in %s'%(n_obj, bin_key))
        num_density = get_surface_number_density(n_objs=n_obj, area_in_sr=area_in_sr)
        N[:, i, i] = 1./num_density

    ###############################################################################################
    # plot all cls
    if plot_cls:
        plt.clf()
        nrow, ncol = 2, 2
        fig, axes = plt.subplots(nrow, ncol)
        # plot cross spectra
        for i, bin1_key in enumerate(all_bins):
            # plot auto spectrum and noise
            p = axes[0, 0].plot(ell, S[:, i, i], 'o-', label=r'$%s$'%(bin1_key))
            axes[1, 0].plot(ell, N[:, i, i], 'o:', color= p[0].get_color(), label='$%s$'%(bin1_key))
            # plot cross spectrum
            for j, bin2_key in enumerate(all_bins):
                if i!=j: # cross spectrum
                    if j>i:
                        p = axes[0, 1].plot(ell, S[:, i, j], 'o-', label=r'$%s$ - $%s$'%(bin1_key, bin2_key))
                        axes[1, 1].plot(ell, N[:, i, j], 'o:', color=p[0].get_color(), label='$%s$ - $%s$'%(bin1_key, bin2_key))
                    else:
                        # check to ensure that we really dont need to plot j, i entry.
                        if (abs(S[:, i, j]-S[:, j, i])>1e-5).any():
                            print('Somethings wrong. Signal matrix isnt symmetric: S[:, i, j]-S[:, j, i]:\n%s'%(S[:, i, j]-S[:, j, i]))
                        if (abs(N[:, i, j]-N[:, j, i])>1e-5).any():
                            print('Somethings wrong. Noise matrix isnt symmetric: N[:, i, j]-N[:, j, i]:\n%s'%(N[:, i, j]-N[:, j, i]))
        for row in range(nrow):
            for col in range(ncol):
                axes[1, col].set_xlabel('$\ell$', fontsize=fontsize)
                axes[row, col].ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
                axes[row, col].set_xscale('log')
                axes[0, col].set_yscale('log')
                axes[row, col].tick_params(axis='x', labelsize=fontsize-2)
                axes[row, col].tick_params(axis='y', labelsize=fontsize-2)
        axes[0, 0].set_ylabel('Auto C$_\ell$', fontsize=fontsize)
        axes[1, 0].set_ylabel('Auto Noise', fontsize=fontsize)
        axes[0, 1].set_ylabel('Cross C$_\ell$', fontsize=fontsize)
        axes[1, 1].set_ylabel('Cross Noise', fontsize=fontsize)
        axes[0, 0].legend(bbox_to_anchor=(-0.15,1), fontsize=fontsize-4)
        axes[0, 1].legend(bbox_to_anchor=(1,1), fontsize=fontsize-4)
        plt.suptitle(filetag, fontsize=fontsize)
        plt.gcf().set_size_inches(20, 10)
        if save_plots:
            filename = '%s_dndz_%sbins_spectra.png'%(filetag, len(z_phot)-1)
            plt.savefig('%s/%s'%(outDir, filename), format='png', bbox_inches='tight')
            print('\n## Saved plot: %s\n'%filename)
        if dont_show_plots:
            plt.close('all')
        else:
            plt.show()

    ###############################################################################################
    # (2*ell+1) Tr(S_ell . C_ell^-1 . S_ell . C_ell^-1)
    C = S+N
    sn_sq = 0
    for j in range(len(ell)):
        inv = np.linalg.inv(C[j, :, :])
        mat = np.dot(S[j, :, :], inv)
        mat = np.dot(inv, mat)
        mat = np.dot(S[j, :, :], mat)
        sn_sq += (2*ell[j]+1)*np.trace(mat)

    fsky = area_in_sr/(4*np.pi)  # total sky area: 4pi Sr
    print('\n## fsky: %.2e\n'%fsky)
    return np.sqrt((fsky/2.)*sn_sq)