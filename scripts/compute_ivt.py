'''Compute Integrated Vapor Transport (IVT) from u- and v- vapor flux components.

IVT = sqrt(uflux^2 + vflux^2), in kg/m/s.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-03-28 12:09:59.
'''

#--------------Globals------------------------------------------
#-----------uflux----------------------
UFLUX_FILE='/home/guangzhi/datasets/erai_qflux/uflux_m1-60_6_2007_cln-cea-proj.nc'
UFLUX_VARID='uflux'

#-----------vflux----------------------
VFLUX_FILE='/home/guangzhi/datasets/erai_qflux/vflux_m1-60_6_2007_cln-cea-proj.nc'
VFLUX_VARID='vflux'

OUTPUTFILE='/home/guangzhi/datasets/quicksave2/THR/ivt_m1-60_6_2007_crop2.nc';



#--------Import modules-------------------------
import numpy as np
from ipart.utils import funcs2 as funcs
from ipart.utils import plot2 as plot



#-------------Main---------------------------------
if __name__=='__main__':

    #-----------Read in data----------------------
    uflux=funcs.readNC(UFLUX_FILE, UFLUX_VARID)
    vflux=funcs.readNC(VFLUX_FILE, VFLUX_VARID)

    ivtdata=np.ma.sqrt(uflux.data**2+vflux.data**2)
    ivt=funcs.NCVAR(ivtdata, 'ivt', uflux.axislist, {'name': 'ivt',
        'long_name': 'integrated vapor transport (IVT)',
        'standard_name': 'integrated vapor transport (IVT)',
        'title': 'integrated vapor transport (IVT)',
        'units': getattr(uflux, 'units', '')})

    #--------Save------------------------------------
    print('\n### <compute_ivt>: Saving output to:\n',OUTPUTFILE)
    funcs.saveNC(OUTPUTFILE, ivt, 'w')

    #-------------------Plot------------------------
    import matplotlib.pyplot as plt
    figure=plt.figure(figsize=(7,10),dpi=100)
    idx=40
    time_str=str(uflux.getTime()[idx])

    plot_vars=[uflux.data[idx], vflux.data[idx], ivt.data[idx]]
    titles=['U', 'V', 'IVT']

    for ii, vii in enumerate(plot_vars):
        axii=figure.add_subplot(3,1,ii+1)
        iso=plot.Isofill(vii, 10, 1, 2)
        plot.plot2(vii, iso, axii,
                title='%s, %s' %(str(time_str), titles[ii]),
                xarray=uflux.getLongitude(),
                yarray=uflux.getLatitude(),
                legend='local',
                projection='cyl',
                fix_aspect=False)

    figure.show()


