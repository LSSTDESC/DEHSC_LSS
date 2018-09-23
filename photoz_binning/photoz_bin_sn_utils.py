import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt

__all__ = ['get_dn_dz', 'get_bins_ends', 'get_bin_edges',
          'get_surface_number_density', 'calc_sn']

# --------------------------------------------------------------------------------------------------------
def get_dn_dz(z_bins, hsc_z_phot, hsc_ids, matched_pdf_ids, matched_pdfs):
    # to get dn/dz, we first find all the galaxies whose z_phot are in a given bin
    # then dn/dz is the sum of the pdfs of these galaxies
    print('Running get_dn_dz ... ')
    dndz = {}
    for i in range(len(z_bins)-1):
        zmin, zmax = z_bins[i], z_bins[i+1]
        bin_key = '%s<z<%s'%(zmin, zmax)
        dndz[bin_key] = []
        ind = np.where((hsc_z_phot >= zmin) & (hsc_z_phot < zmax))[0]
        print('Considering pdfs of %s objects for %s'%(len(ind), bin_key))

        for j, ID in enumerate(hsc_ids[ind]):
            row = np.where(matched_pdf_ids==ID)[0][0]
            if j==0:
                dndz[bin_key] = matched_pdfs[row]
            else:
                dndz[bin_key] += matched_pdfs[row]
    return dndz

# --------------------------------------------------------------------------------------------------------
def get_bins_ends(z_bin_array):
    # get the labels and bin edges based on the input z-bin array
    bin_keys, bin_ends = [], {}
    for i in range(len(z_bin_array)-1):
        zmin, zmax = z_bin_array[i], z_bin_array[i+1]
        bin_key = '%s<z<%s'%(zmin, zmax)
        bin_keys.append(bin_key)
        bin_ends[bin_key] = [zmin, zmax]
    return np.array(bin_keys), bin_ends

# --------------------------------------------------------------------------------------------------------
def get_bin_edges(nbin, hsc_z_phot, z_bins):
    # decide on bins edges: we look for bin edges that have roughly the same number of galaxies
    zmin, zmax = 0.15, 1.5   # these are fixed since this is the range where z_phot is reliable
    obj_thres = 3000    # the bins would have galaxies within this threshold of each other
    if nbin==0:
        print('nbin must be greater than one.')
        return
    else:
        n_obj = len(np.where((hsc_z_phot > zmin) & (hsc_z_phot < zmax))[0])
        if nbin==1:
            print('%s objects in %s<z<%s'%(n_obj, zmin, zmax))
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
                good_edges = [f for f in z_bins if f>bin_ends[i]]

                while abs(n_obj-wanted_n_obj_in_bin)>obj_thres and j<len(good_edges):
                    z_edge = good_edges[j]
                    n_obj = len(np.where((hsc_z_phot > bin_ends[i]) & (hsc_z_phot < z_edge))[0])
                    #print('## %s objects in %s<z<%s'%(n_obj, bin_ends[i], z_edge))
                    j += 1
                if abs(n_obj-wanted_n_obj_in_bin)>obj_thres:
                    print('Something is wrong.')

                if (z_edge==bin_ends[-1]):
                    print('Need to change the threshold on the number of galaxies: %s currently'%(obj_thres))

                bin_ends[i+1] = float('%.2f'%z_edge)

            print('\nFinal:')
            for i in range(nbin):
                n_obj = len(np.where((hsc_z_phot > bin_ends[i]) & (hsc_z_phot < bin_ends[i+1]))[0])
                print('%s objects in %s<z<%s'%(n_obj, bin_ends[i], bin_ends[i+1]))

        return bin_ends
# --------------------------------------------------------------------------------------------------------
def get_surface_number_density(n_objs, area_in_sr):  # in 1/Sr
    return n_objs/area_in_sr

# --------------------------------------------------------------------------------------------------------
def calc_sn(z_phot, z_bins, hsc_z_phot, hsc_ids, matched_pdf_ids, matched_pdfs, n_z, ell, area_in_sr, plot_cls=True):
    print('-----------------------------------------\nRunning calc_sn ... ')

    z_all = [min(z_bins)] + list(z_phot) + [max(z_bins)]   # cover the entire range, even outside the target bins

    target_bins, target_bin_ends = get_bins_ends(z_bin_array=z_phot)
    all_bins, all_bin_ends = get_bins_ends(z_bin_array=z_all)

    print('Target bins: %s'%(target_bins))
    print('All bins: %s\n'%(all_bins))

    dn_dz = get_dn_dz(z_bins=z_all, hsc_z_phot=hsc_z_phot, hsc_ids=hsc_ids,
                      matched_pdf_ids=matched_pdf_ids, matched_pdfs=matched_pdfs)
    plt.clf()
    for key in dn_dz:
        if key in target_bin_ends:
            plt.plot(z_bins, dn_dz[key], lw=4, label='** dn_dz: %s'%key)
            for z_edge in target_bin_ends[key]:
                plt.plot([z_edge, z_edge], [0,7000], 'k-.', alpha=0.7)
        else:
            plt.plot(z_bins, dn_dz[key], lw=4, label='dn_dz: %s'%key)
    plt.plot(z_bins, n_z, 'k:', label='n_z', alpha=0.7)
    plt.legend()
    plt.xlabel('z')
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
        for i, bin1_key in enumerate(all_bins):
            for j, bin2_key in enumerate(all_bins):
                if i!=j:
                    plt.plot(ell, S[:, i, j], 'o-', label='%s-%s'%(bin1_key, bin2_key))
        plt.xlabel('$\ell$')
        plt.ylabel('Cross C$_\ell$')
        plt.gca().set_xscale('log')
        plt.gca().ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
        plt.legend(bbox_to_anchor=(1,.5))
        plt.show()

        # plot all cls
        for i, bin1_key in enumerate(all_bins):
            for j, bin2_key in enumerate(all_bins):
                if i==j:
                    plt.plot(ell, S[:, i, j], 'o-', label='%s-%s'%(bin1_key, bin2_key))
        plt.xlabel('$\ell$')
        plt.ylabel('Auto C$_\ell$')
        plt.gca().ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
        plt.gca().set_xscale('log')
        plt.legend(bbox_to_anchor=(1,.5))
        plt.show()

        for i, bin1_key in enumerate(all_bins):
            for j, bin2_key in enumerate(all_bins):
                plt.plot(ell, N[:, i, j], 'o:', label='%s-%s'%(bin1_key, bin2_key))
        plt.xlabel('$\ell$')
        plt.ylabel('Noise C$_\ell$')
        plt.gca().ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
        plt.gca().set_xscale('log')
        plt.legend(bbox_to_anchor=(1,.5))
        plt.show()

    ###############################################################################################
    # (2*ell+1) Tr(S_ell . C_ell^-1 . S_ell . C_ell^-1)
    C = S+N
    sn = 0
    for j in range(len(ell)):
        inv = np.linalg.inv(C[j, :, :])
        mat = np.dot(S[j, :, :], inv)
        mat = np.dot(inv, mat)
        mat = np.dot(S[j, :, :], mat)
        sn += (2*ell[j]+1)*np.trace(mat)

    return sn