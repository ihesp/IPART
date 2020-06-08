'''Compute Integrated Vapor Transport (IVT) from u- and v- vapor flux components.

IVT = sqrt(uflux^2 + vflux^2), in kg/m/s.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-03-28 12:09:59.
'''

#--------------Globals------------------------------------------
#-----------uflux----------------------
UFLUX_FILE='/home/guangzhi/datasets/erai/ERAI_AR_THR/uflux_m1-60_6_1984_crop.nc'
UFLUX_VARID='uflux'

#-----------vflux----------------------
VFLUX_FILE='/home/guangzhi/datasets/erai/ERAI_AR_THR/vflux_m1-60_6_1984_crop.nc'
VFLUX_VARID='vflux'

OUTPUTFILE='/home/guangzhi/datasets/erai/ERAI_AR_THR/ivt_m1-60_6_1984_crop2.nc';



#--------Import modules-------------------------
import cdms2 as cdms
import MV2 as MV
from utils import funcs, plot



#-------------Main---------------------------------
if __name__=='__main__':

    #-----------Read in data----------------------
    uflux=funcs.readVar(UFLUX_FILE, UFLUX_VARID)
    vflux=funcs.readVar(VFLUX_FILE, VFLUX_VARID)

    ivt=MV.sqrt(uflux*uflux+vflux*vflux)
    ivt.id='ivt'
    ivt.long_name='integrated vapor transport (IVT)'
    ivt.standard_name=ivt.long_name
    ivt.title=ivt.long_name
    ivt.units=getattr(uflux, 'units', '')

    #--------Save------------------------------------
    print('\n### <compute_ivt>: Saving output to:\n',OUTPUTFILE)
    #fout=cdms.open(OUTPUTFILE,'w')
    #fout.write(ivt,typecode='f')
    #fout.close()

    #-------------------Plot------------------------
    import matplotlib.pyplot as plt
    figure=plt.figure(figsize=(7,10),dpi=100)
    idx=40
    time_str=uflux.getTime().asComponentTime()[idx]

    plot_vars=[uflux[idx], vflux[idx], ivt[idx]]
    titles=['U', 'V', 'IVT']

    for ii, vii in enumerate(plot_vars):
        axii=figure.add_subplot(3,1,ii+1)
        iso=plot.Isofill(vii, 10, 1, 2)
        plot.plot2(vii, iso, axii,
                title='%s, %s' %(str(time_str), titles[ii]),
                legend='local',
                projection='cyl',
                fix_aspect=False)

    figure.show()


