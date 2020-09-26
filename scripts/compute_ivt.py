'''Compute Integrated Vapor Transport (IVT) from u- and v- vapor flux components.

u-vapor flux component's standard_name:
    "eastward_atmosphere_water_transport_across_unit_distance".
v-vapor flux component's standard_name:
    "northward_atmosphere_water_transport_across_unit_distance".

IVT = sqrt(uflux^2 + vflux^2), in kg/m/s.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 10:08:31.
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
from ipart.utils import funcs
from ipart.utils import plot



#-------------Main---------------------------------
if __name__=='__main__':

    #-----------Read in data----------------------
    ufluxNV=funcs.readNC(UFLUX_FILE, UFLUX_VARID)
    vfluxNV=funcs.readNC(VFLUX_FILE, VFLUX_VARID)

    ivtdata=np.ma.sqrt(ufluxNV.data**2+vfluxNV.data**2)
    ivtNV=funcs.NCVAR(ivtdata, 'ivt', ufluxNV.axislist, {'name': 'ivt',
        'long_name': 'integrated vapor transport (IVT)',
        'standard_name': 'integrated vapor transport (IVT)',
        'title': 'integrated vapor transport (IVT)',
        'units': getattr(ufluxNV, 'units', '')})

    #--------Save------------------------------------
    print('\n### <compute_ivt>: Saving output to:\n',OUTPUTFILE)
    funcs.saveNC(OUTPUTFILE, ivtNV, 'w')

    #-------------------Plot------------------------
    import matplotlib.pyplot as plt
    figure=plt.figure(figsize=(7,10),dpi=100)
    idx=40
    time_str=str(ufluxNV.getTime()[idx])

    plot_vars=[ufluxNV.data[idx], vfluxNV.data[idx], ivtNV.data[idx]]
    titles=['U', 'V', 'IVT']

    for ii, vii in enumerate(plot_vars):
        axii=figure.add_subplot(3,1,ii+1)
        iso=plot.Isofill(vii, 10, 1, 2)
        plot.plot2(vii, iso, axii,
                title='%s, %s' %(str(time_str), titles[ii]),
                xarray=ufluxNV.getLongitude(),
                yarray=ufluxNV.getLatitude(),
                legend='local',
                fix_aspect=False)

    figure.show()


