import numpy as np
import os
import sys

def submit_job(fname_sql,output_format,do_preview=False) :
    command="python hscReleaseQuery.py --user=damonge@local"
    if do_preview :
        command+=" -p"
    command+=" -f "+output_format
    command+=" "+fname_sql
    print(command)
    os.system(command)

def write_boxsearch(tablename,dec_range,ra_range,fname_out,output_format="fits",do_box=True,submit=False,strict_cuts=False) :
    filters=['g','r','i','z','y']
    stout="-- Run box, "+fname_out+"\n"

    def add_filters(name,behind=True) :
        sthere=""
        for fi in filters :
            if behind :
                sthere+="       ,"+name+fi+"\n"
            else :
                sthere+="       ,"+fi+name+"\n"
        return sthere

    stout+="SELECT object_id\n"
    stout+="       ,ra\n"
    stout+="       ,dec\n"
    stout+="       ,tract\n"
    stout+="       ,patch\n"
    stout+=add_filters("merge_peak_")
    stout+=add_filters("countinputs",behind=False)
    stout+="       ,iflags_pixel_bright_object_center\n"
    stout+="       ,iflags_pixel_bright_object_any\n" 
    stout+="       ,iclassification_extendedness\n"
    #Dust extinction
    stout+=add_filters("a_")
    #Psf fluxes and magnitudes
    stout+=add_filters("flux_psf",behind=False)
    stout+=add_filters("flux_psf_err",behind=False)
    stout+=add_filters("flux_psf_flags",behind=False)
    stout+=add_filters("mag_psf",behind=False)
    stout+=add_filters("mag_psf_err",behind=False)
    #Aperture fluxes and magnitudes
    stout+=add_filters("flux_aperture10",behind=False)
    stout+=add_filters("flux_aperture10_err",behind=False)
    stout+=add_filters("flux_aperture_flags",behind=False)
    stout+=add_filters("mag_aperture10",behind=False)
    stout+=add_filters("mag_aperture10_err",behind=False)
    #Cmodel fluxes and magnitudes
    stout+=add_filters("cmodel_flux",behind=False)
    stout+=add_filters("cmodel_flux_err",behind=False)
    stout+=add_filters("cmodel_flux_flags",behind=False)
    stout+=add_filters("cmodel_mag",behind=False)
    stout+=add_filters("cmodel_mag_err",behind=False)
    stout+="FROM\n"
    stout+="       "+tablename+"\n"
    stout+="WHERE\n"
    stout+="       detect_is_primary=True and\n"
    stout+="       icmodel_flags_badcentroid=False and\n"
    stout+="       icentroid_sdss_flags=False and\n"
    stout+="       iflags_pixel_edge=False and\n"
    stout+="       iflags_pixel_interpolated_center=False and\n"
    stout+="       iflags_pixel_saturated_center=False and\n"
    stout+="       iflags_pixel_cr_center=False and\n"
    stout+="       iflags_pixel_bad=False and\n"
    stout+="       iflags_pixel_suspect_center=False and\n"
    stout+="       iflags_pixel_clipped_any=False\n"
    if strict_cuts :
        stout+=" and\n"
        stout+="       gcentroid_sdss_flags=False and\n"
        stout+="       rcentroid_sdss_flags=False and\n"
        #stout+="       icentroid_sdss_flags=False and\n"
        stout+="       zcentroid_sdss_flags=False and\n"
        stout+="       ycentroid_sdss_flags=False and\n"
        stout+="       gcmodel_flux_flags=False and\n"
        stout+="       rcmodel_flux_flags=False and\n"
        stout+="       icmodel_flux_flags=False and\n"
        stout+="       zcmodel_flux_flags=False and\n"
        stout+="       ycmodel_flux_flags=False and\n"
        stout+="       gflux_psf_flags=False and\n"
        stout+="       rflux_psf_flags=False and\n"
        stout+="       iflux_psf_flags=False and\n"
        stout+="       zflux_psf_flags=False and\n"
        stout+="       yflux_psf_flags=False and\n"
        stout+="       gflags_pixel_edge=False and\n"
        stout+="       rflags_pixel_edge=False and\n"
        #stout+="       iflags_pixel_edge=False and\n"
        stout+="       zflags_pixel_edge=False and\n"
        stout+="       yflags_pixel_edge=False and\n"
        stout+="       gflags_pixel_interpolated_center=False and\n"
        stout+="       rflags_pixel_interpolated_center=False and\n"
        #stout+="       iflags_pixel_interpolated_center=False and\n"
        stout+="       zflags_pixel_interpolated_center=False and\n"
        stout+="       yflags_pixel_interpolated_center=False and\n"
        stout+="       gflags_pixel_saturated_center=False and\n"
        stout+="       rflags_pixel_saturated_center=False and\n"
        #stout+="       iflags_pixel_saturated_center=False and\n"
        stout+="       zflags_pixel_saturated_center=False and\n"
        stout+="       yflags_pixel_saturated_center=False and\n"
        stout+="       gflags_pixel_cr_center=False and\n"
        stout+="       rflags_pixel_cr_center=False and\n"
        #stout+="       iflags_pixel_cr_center=False and\n"
        stout+="       zflags_pixel_cr_center=False and\n"
        stout+="       yflags_pixel_cr_center=False and\n"
        stout+="       gflags_pixel_bad=False and\n"
        stout+="       rflags_pixel_bad=False and\n"
        #stout+="       iflags_pixel_bad=False and\n"
        stout+="       zflags_pixel_bad=False and\n"
        stout+="       yflags_pixel_bad=False\n"
    if do_box :
        stout+=" and\n"
        stout+="       boxSearch(coord,%.3lf,%.3lf,%.3lf,%.3lf)\n"%(ra_range[0],ra_range[1],dec_range[0],dec_range[1])
    else :
        stout+=" \n"
    stout+=";\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

def add_photoz(mt) :
    sto=""
    sto+="       ,p"+mt+".photoz_mean as pz_mean_"+mt+"\n"
    sto+="       ,p"+mt+".photoz_mode as pz_mode_"+mt+"\n"
    sto+="       ,p"+mt+".photoz_best as pz_best_"+mt+"\n"
    sto+="       ,p"+mt+".photoz_mc as pz_mc_"+mt+"\n"
    return sto

def write_boxsearch_random(tablename,dec_range,ra_range,fname_out,output_format="fits",do_box=True,submit=False) :
    filters=['g','r','i','z','y']
    stout="-- Run box, "+fname_out+"\n"

    def add_filters(name,behind=True) :
        sthere=""
        for fi in filters :
            if behind :
                sthere+="       ,"+name+fi+"\n"
            else :
                sthere+="       ,"+fi+name+"\n"
        return sthere

    stout+="SELECT object_id\n"
    stout+="       ,ra\n"
    stout+="       ,dec\n"
    stout+="       ,tract\n"
    stout+="       ,patch\n"
    stout+="       ,idetect_is_primary\n"
    stout+="       ,shape_detradius(array[ishape_sdss_psf_11, ishape_sdss_psf_22, ishape_sdss_psf_12])*2.0*sqrt(2.0)*0.17 as ipsf_size\n"
    stout+=add_filters("countinputs",behind=False)
    stout+="       ,iflags_pixel_bright_object_center\n"
    stout+="       ,iflags_pixel_bright_object_any\n" 
    stout+="       ,iadjust_density\n"
    #Psf fluxes and magnitudes
    stout+=add_filters("sky_mean",behind=False)
    stout+=add_filters("sky_std",behind=False)
    stout+=add_filters("pix_variance",behind=False)
    stout+="FROM\n"
    stout+="       "+tablename+"\n"
    stout+="WHERE\n"
    stout+="       iflags_pixel_edge=False and\n"
    stout+="       iflags_pixel_interpolated_center=False and\n"
    stout+="       iflags_pixel_saturated_center=False and\n"
    stout+="       iflags_pixel_cr_center=False and\n"
    stout+="       iflags_pixel_bad=False and\n"
    stout+="       iflags_pixel_suspect_center=False and\n"
    stout+="       ra between %.3lf and %.3lf and dec between %.3lf and %.3lf\n"%(ra_range[0],ra_range[1],dec_range[0],dec_range[1])
    stout+=";\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

def write_frames(tablename,fname_out,output_format='fits',submit=False) :
    stout="-- Run metadata\n"
    stout+="SELECT *\n"
    stout+="FROM "+tablename+".frame\n"
    stout+=";\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

def write_fieldsearch(tablename,fieldname,fname_out,output_format="fits",do_field=True,submit=False,ra_range=None,do_photoz=False,strict_cuts=False) :
    filters=['g','r','i','z','y']
    stout="-- Run field, "+fname_out+"\n"

    def add_filters(name,behind=True) :
        sthere=""
        for fi in filters :
            if behind :
                sthere+="       ,forced."+name+fi+"\n"
            else :
                sthere+="       ,forced."+fi+name+"\n"
        return sthere

    stout+="SELECT object_id\n"
    stout+="       ,forced.ra as ra\n"
    stout+="       ,forced.dec as dec\n"
    stout+="       ,forced.tract as tract\n"
    stout+="       ,forced.patch as patch\n"
    stout+=add_filters("merge_peak_")
    stout+=add_filters("countinputs",behind=False)
    stout+="       ,forced.iflags_pixel_bright_object_center\n"
    stout+="       ,forced.iflags_pixel_bright_object_any\n" 
    stout+="       ,forced.iclassification_extendedness\n"
    stout+="       ,meas.iblendedness_abs_flux as iblendedness_abs_flux\n"
    #Dust extinction
    stout+=add_filters("a_")
    #Psf fluxes and magnitudes
    stout+=add_filters("flux_psf",behind=False)
    stout+=add_filters("flux_psf_err",behind=False)
    stout+=add_filters("flux_psf_flags",behind=False)
    stout+=add_filters("mag_psf",behind=False)
    stout+=add_filters("mag_psf_err",behind=False)
    #Aperture fluxes and magnitudes
    stout+=add_filters("flux_aperture10",behind=False)
    stout+=add_filters("flux_aperture10_err",behind=False)
    stout+=add_filters("flux_aperture_flags",behind=False)
    stout+=add_filters("mag_aperture10",behind=False)
    stout+=add_filters("mag_aperture10_err",behind=False)
    #Cmodel fluxes and magnitudes
    stout+=add_filters("cmodel_flux",behind=False)
    stout+=add_filters("cmodel_flux_err",behind=False)
    stout+=add_filters("cmodel_flux_flags",behind=False)
    stout+=add_filters("cmodel_mag",behind=False)
    stout+=add_filters("cmodel_mag_err",behind=False)
    if do_photoz :
        stout+=add_photoz("eab")
        stout+=add_photoz("frz")
        stout+=add_photoz("nnz")
    stout+="FROM\n"
    stout+="       "+tablename+".forced as forced\n"
    if do_photoz :
        stout+="       LEFT JOIN "+tablename+".meas meas USING (object_id)\n"
        stout+="       LEFT JOIN "+tablename+".photoz_ephor_ab peab USING (object_id)\n"
        stout+="       LEFT JOIN "+tablename+".photoz_frankenz pfrz USING (object_id)\n"
        stout+="       LEFT JOIN "+tablename+".photoz_nnpz pnnz USING (object_id)\n"
    stout+="WHERE\n"
    stout+="       forced.detect_is_primary=True and\n"
    stout+="       forced.icmodel_flags_badcentroid=False and\n"
    stout+="       forced.icentroid_sdss_flags=False and\n"
    stout+="       forced.iflags_pixel_edge=False and\n"
    stout+="       forced.iflags_pixel_interpolated_center=False and\n"
    stout+="       forced.iflags_pixel_saturated_center=False and\n"
    stout+="       forced.iflags_pixel_cr_center=False and\n"
    stout+="       forced.iflags_pixel_bad=False and\n"
    stout+="       forced.iflags_pixel_suspect_center=False and\n"
    stout+="       forced.iflags_pixel_clipped_any=False and\n"
    stout+="       meas.ideblend_skipped=False"
    if strict_cuts :
        stout+=" and\n"
        stout+="       forced.gcentroid_sdss_flags=False and\n"
        stout+="       forced.rcentroid_sdss_flags=False and\n"
        #stout+="       forced.icentroid_sdss_flags=False and\n"
        stout+="       forced.zcentroid_sdss_flags=False and\n"
        stout+="       forced.ycentroid_sdss_flags=False and\n"
        stout+="       forced.gcmodel_flux_flags=False and\n"
        stout+="       forced.rcmodel_flux_flags=False and\n"
        stout+="       forced.icmodel_flux_flags=False and\n"
        stout+="       forced.zcmodel_flux_flags=False and\n"
        stout+="       forced.ycmodel_flux_flags=False and\n"
        stout+="       forced.gflux_psf_flags=False and\n"
        stout+="       forced.rflux_psf_flags=False and\n"
        stout+="       forced.iflux_psf_flags=False and\n"
        stout+="       forced.zflux_psf_flags=False and\n"
        stout+="       forced.yflux_psf_flags=False and\n"
        stout+="       forced.gflags_pixel_edge=False and\n"
        stout+="       forced.rflags_pixel_edge=False and\n"
        #stout+="       forced.iflags_pixel_edge=False and\n"
        stout+="       forced.zflags_pixel_edge=False and\n"
        stout+="       forced.yflags_pixel_edge=False and\n"
        stout+="       forced.gflags_pixel_interpolated_center=False and\n"
        stout+="       forced.rflags_pixel_interpolated_center=False and\n"
        #stout+="       forced.iflags_pixel_interpolated_center=False and\n"
        stout+="       forced.zflags_pixel_interpolated_center=False and\n"
        stout+="       forced.yflags_pixel_interpolated_center=False and\n"
        stout+="       forced.gflags_pixel_saturated_center=False and\n"
        stout+="       forced.rflags_pixel_saturated_center=False and\n"
        #stout+="       forced.iflags_pixel_saturated_center=False and\n"
        stout+="       forced.zflags_pixel_saturated_center=False and\n"
        stout+="       forced.yflags_pixel_saturated_center=False and\n"
        stout+="       forced.gflags_pixel_cr_center=False and\n"
        stout+="       forced.rflags_pixel_cr_center=False and\n"
        #stout+="       forced.iflags_pixel_cr_center=False and\n"
        stout+="       forced.zflags_pixel_cr_center=False and\n"
        stout+="       forced.yflags_pixel_cr_center=False and\n"
        stout+="       forced.gflags_pixel_bad=False and\n"
        stout+="       forced.rflags_pixel_bad=False and\n"
        #stout+="       forced.iflags_pixel_bad=False and\n"
        stout+="       forced.zflags_pixel_bad=False and\n"
        stout+="       forced.yflags_pixel_bad=False\n"
# and\n"
#    stout+="       iclassification_extendedness!=0"
    if do_field :
        stout+=" and\n"
        stout+="       "+tablename+".search_"+fieldname+"(object_id)"
    if ra_range is not None :
        stout+=" and\n"
        stout+="       forced.ra>=%.3lf and\n"%(ra_range[0])
        stout+="       forced.ra<%.3lf"%(ra_range[1])
    stout+="\n;\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

def write_countsget(tablename,dec_range,ra_range,fname_out,output_format='csv',do_box=True,submit=False,is_random=False) :
    stout="-- Run count, "+fname_out+"\n"
    stout+="SELECT count(*)\n"
    stout+="FROM\n"
    stout+="       "+tablename+"\n"
    if do_box or (not is_random) :
        stout+="WHERE\n"
        if not is_random :
            stout+="       detect_is_primary=True and\n"
            stout+="       icmodel_flags_badcentroid=False and\n"
            stout+="       icentroid_sdss_flags=False and\n"
            stout+="       iflags_pixel_edge=False and\n"
            stout+="       iflags_pixel_interpolated_center=False and\n"
            stout+="       iflags_pixel_saturated_center=False and\n"
            stout+="       iflags_pixel_cr_center=False and\n"
            stout+="       iflags_pixel_bad=False and\n"
            stout+="       iflags_pixel_suspect_center=False and\n"
            stout+="       iflags_pixel_clipped_any=False and\n"
        if do_box :
            if is_random :
                stout+="ra BETWEEN %.3lf AND %.3lf and dec between %.3lf and %.3lf\n"%(ra_range[0],ra_range[1],dec_range[0],dec_range[1])
            else :
                stout+="       boxSearch(coord,%.3lf,%.3lf,%.3lf,%.3lf)\n"%(ra_range[0],ra_range[1],dec_range[0],dec_range[1])
    stout+=";\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

def write_getspecz(tablename,fname_out,output_format='fits',submit=False) :
    stout="-- Run spec, "+fname_out+"\n"
    stout+="SELECT *\n"
    stout+="FROM\n"
    stout+="       "+tablename+"\n"
    stout+=";\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

#Per-frame metadata
write_frames("pdr1_wide","frames_wide.sql",submit=True)
write_frames("pdr1_deep","frames_deep.sql",submit=True)
write_frames("pdr1_udeep","frames_udeep.sql",submit=True)

#WIDE fields
for fld in ['aegis','gama09h','gama15h','hectomap','wide12h','xmm_lss'] :
    write_fieldsearch("pdr1_wide",fld,"field_wide_"+fld+"_pz.sql",do_field=True,submit=True,do_photoz=True)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h1_pz.sql",do_field=True,submit=True,ra_range=[330.,336.],do_photoz=True)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h2_pz.sql",do_field=True,submit=True,ra_range=[336.,342.],do_photoz=True)

#WIDE fields, strict cuts
for fld in ['gama09h','gama15h','hectomap','wide12h','xmm_lss'] :
    write_fieldsearch("pdr1_wide",fld,"field_wide_"+fld+"_pz_strict.sql",do_field=True,submit=True,do_photoz=True,strict_cuts=True)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h1_pz_strict.sql",do_field=True,submit=True,ra_range=[330.,336.],do_photoz=True,strict_cuts=True)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h2_pz_strict.sql",do_field=True,submit=True,ra_range=[336.,342.],do_photoz=True,strict_cuts=True)

#DEEP fields
for fld in ['cosmos','deep2_3','elais_n1','xmm_lss'] :
    write_fieldsearch("pdr1_deep",fld,"field_deep_"+fld+".sql",do_field=True,submit=True,do_photoz=True)

#UDEEP fields
for fld in ['cosmos','sxds'] :
    write_fieldsearch("pdr1_udeep",fld,"field_udeep_"+fld+".sql",do_field=True,submit=True,do_photoz=True)
exit(1)
#WIDE-depth COSMOS
for see in ['best','median','worst'] :
    write_fieldsearch("pdr1_cosmos_widedepth_"+see,"none","field_cosmo_wide_"+see+".sql",do_field=False,submit=True,do_photoz=False)

#Frames table


exit(1)

#WIDE fields
for fld in ['aegis','gama09h','gama15h','hectomap','wide12h','xmm_lss'] :
    write_fieldsearch("pdr1_wide",fld,"field_wide_"+fld+".sql",do_field=True,submit=True)
#VVDS is too large for a single query
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h1.sql",do_field=True,submit=True,ra_range=[330.,336.])
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h2.sql",do_field=True,submit=True,ra_range=[336.,342.])

#WIDE-depth COSMOS
for see in ['best','median','worst'] :
    write_fieldsearch("pdr1_cosmos_widedepth_"+see,"none","field_cosmo_wide_"+see+".sql",do_field=False,submit=True)

#DEEP fields
for fld in ['cosmos','deep2_3','elais_n1','xmm_lss'] :
    write_fieldsearch("pdr1_deep",fld,"field_deep_"+fld+".sql",do_field=True,submit=True)

#UDEEP fields
for fld in ['cosmos','sxds'] :
    write_fieldsearch("pdr1_udeep",fld,"field_udeep_"+fld+".sql",do_field=True,submit=True)
