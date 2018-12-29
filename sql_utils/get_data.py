from pdr1_queries import write_frames, write_fieldsearch

#Per-frame metadata
write_frames("pdr1_wide","frames_wide.sql",submit=True)
write_frames("pdr1_deep","frames_deep.sql",submit=True)
write_frames("pdr1_udeep","frames_udeep.sql",submit=True)

#WIDE fields
for fld in ['aegis','gama09h','gama15h','hectomap','wide12h','xmm_lss'] :
    write_fieldsearch("pdr1_wide",fld,"field_wide_"+fld+"_pz_strict.sql",do_field=True,
                      submit=True,do_photoz=True,strict_cuts=True)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h1_pz_strict.sql",do_field=True,
                  submit=True,ra_range=[330.,336.],do_photoz=True,strict_cuts=True,part=1)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h2_pz_strict.sql",do_field=True,
                  submit=True,ra_range=[336.,342.],do_photoz=True,strict_cuts=True,part=2)

#DEEP fields
for fld in ['cosmos','deep2_3','elais_n1','xmm_lss'] :
    write_fieldsearch("pdr1_deep",fld,"field_deep_"+fld+"_pz_strict.sql",do_field=True,
                      submit=True,do_photoz=True,strict_cuts=True)

#UDEEP fields
for fld in ['cosmos','sxds'] :
    write_fieldsearch("pdr1_udeep",fld,"field_udeep_"+fld+"_pz_strict.sql",do_field=True,
                      submit=True,do_photoz=True,strict_cuts=True)

#WIDE-depth COSMOS
for see in ['best','median','worst'] :
    write_fieldsearch("pdr1_cosmos_widedepth_"+see,"none","field_cosmo_wide_"+see+".sql",
                      do_field=False,submit=True,do_photoz=False)

#WIDE fields, shear catalog
for fld in ['aegis','gama09h','gama15h','hectomap','wide12h','xmm_lss'] :
    write_fieldsearch("pdr1_wide",fld,"field_wide_"+fld+"_pz.sql",do_field=True,
                      submit=True,do_photoz=True)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h1_pz.sql",do_field=True,submit=True,
                  ra_range=[330.,336.],do_photoz=True,part=1)
write_fieldsearch("pdr1_wide",'vvds',"field_wide_vvds_h2_pz.sql",do_field=True,submit=True,
                  ra_range=[336.,342.],do_photoz=True,part=2)
