import numpy as np
import os
import sys
# This script assumes the existence of an environment variable called HSC_USERNAME
# which contains your username to connect to the HSC database
def submit_job(fname_sql,output_format,do_preview=False) :
    command="python hscReleaseQuery.py -u $HSC_USERNAME"
    if do_preview :
        command+=" -p"
    command+=" -f "+output_format
    command+=" "+fname_sql
    print command
    os.system(command)

def write_boxsearch(tablename,dec_range,ra_range,fname_out,output_format="fits",do_box=True,submit=False) :
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

    stout+="SELECT *\n"
    stout+="FROM\n"
    stout+="       "+tablename+"\n"
    if do_box :
        stout+=" AND\n"
        stout+="       boxSearch(coord,%.3lf,%.3lf,%.3lf,%.3lf)\n"%(ra_range[0],ra_range[1],dec_range[0],dec_range[1])
    else :
        stout+=" \n"
    stout+=";\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)


def write_fieldsearch(tablename,fieldname,fname_out,output_format="fits",do_field=True,submit=False,ra_range=None) :

    stout="-- Run field, "+fname_out+"\n"
    stout+="SELECT *\n"
    stout+="FROM\n"
    stout+="       "+tablename+"\n"
    if do_field :
        stout+=" WHERE\n"
        stout+="       "+"pdr1_wide.search_"+fieldname+"(object_id)"
    if ra_range is not None :
        stout+=" WHERE\n"
        stout+="       ra>=%.3lf AND\n"%(ra_range[0])
        stout+="       ra<%.3lf"%(ra_range[1])
    stout+="\n;\n"

    f=open(fname_out,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_out,output_format)

#WIDE fields
for photoz in ['demp','ephor','ephor_ab','frankenz','mizuki','mlz','nnpz']:
    for fld in ['gama09h']:#['aegis','gama09h','gama15h','hectomap','wide12h','xmm_lss'] :
         write_fieldsearch("pdr1_wide.photoz_"+str(photoz),fld,"field_wide_"+fld+"_"+photoz+".sql",do_field=True,submit=True)
#VVDS is too large for a single query
#    write_fieldsearch("pdr1_wide.photoz_"+str(photoz),'vvds',"field_wide_vvds_h1.sql",do_field=True,submit=True,ra_range=[330.,336.])
#    write_fieldsearch("pdr1_wide.photoz_"+str(photoz),'vvds',"field_wide_vvds_h2"+"_"+photoz+".sql",do_field=True,submit=True,ra_range=[336.,342.])

#WIDE-depth COSMOS
#for photoz in ['demp','ephor','ephor_ab','frankenz','mizuki','mlz','nnpz']:
#    for see in ['best','median','worst'] :
#        write_fieldsearch("pdr1_cosmos_widedepth_"+see,"none","field_cosmo_wide_"+see+"_"+photoz+".sql",do_field=False,submit=True)

#DEEP fields
#for photoz in ['demp','ephor','ephor_ab','frankenz','mizuki','mlz','nnpz']:
#    for fld in ['cosmos','deep2_3','elais_n1','xmm_lss'] :
#        write_fieldsearch("pdr1_deep",fld,"field_deep_"+fld+"_"+photoz+".sql",do_field=True,submit=True)

#UDEEP fields
#for photoz in ['demp','ephor','ephor_ab','frankenz','mizuki','mlz','nnpz']:
#    for fld in ['cosmos','sxds'] :
#        write_fieldsearch("pdr1_udeep",fld,"field_udeep_"+fld+"_"+photoz+".sql",do_field=True,submit=True)
