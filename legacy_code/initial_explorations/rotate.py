import numpy as np

DTOR=np.pi/180
RTOD=180/np.pi

def get_rotation_matrix(ra0,dec0) :
    th0=DTOR*(90-dec0)
    ph0=DTOR*ra0
    ct=np.cos(th0); st=np.sin(th0)
    cp=np.cos(ph0); sp=np.sin(ph0)
    rot=np.array([[sp**2+cp**2*st,sp*cp*(st-1),cp*ct],
                  [cp*sp*(st-1),cp**2+sp**2*st,sp*ct],
                  [-cp*ct,-sp*ct,st]])
    return rot

def get_new_coords(data,ra0,dec0) :
    th=DTOR*(90-data['dec'])
    ph=DTOR*data['ra']
    ct=np.cos(th); st=np.sin(th)
    cp=np.cos(ph); sp=np.sin(ph)
    pos=np.array([st*cp,st*sp,ct])
    rot=get_rotation_matrix(ra0,dec0)
    pos=np.dot(rot,pos)
    dec_r=90-RTOD*np.arccos(pos[2])
    ra_r=RTOD*np.arctan2(pos[1],pos[0])
    return rot,ra_r,dec_r
