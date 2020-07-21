'''Check meta data of netcdf file.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-03-31 10:42:41.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os, sys
from netCDF4 import Dataset
from ipart.utils import funcs2 as funcs




#-------------Main---------------------------------
if __name__=='__main__':

    if len(sys.argv)<3:
        raise Exception("Please provide the path to the netcdf file, and the id of the variable.")

    file_path=os.path.abspath(sys.argv[1])

    if not os.path.exists(file_path):
        raise Exception("Input file not exists.")

    print('\n### <test_data>: Read in file:\n',file_path)

    try:
        fin=Dataset(file_path, 'r')
    except Exception as e:
        raise Exception("Failed to open file.\nException: %s" %str(e))
    else:
        fin.close()

    try:
        #var=fin.variables[sys.argv[2]]
        var=funcs.readNC(file_path, sys.argv[2])
    except:
        raise Exception("Variable not found in file. Please verify the variable id is correct.")
    else:

        # check rank
        try:
            if var.ndim not in [3,4]:
                raise Exception("Data should be in rank 3 or 4, i.e. (time, lat, lon) or (time, level, lat, lon).")
        except:
            raise Exception("Data should be in rank 3 or 4, i.e. (time, lat, lon) or (time, level, lat, lon).")


        # try get time axis
        try:
            timeax=var.getTime()
            print('#'*50)
            print('Variable time axis:')
            print('#'*50)
            idx=funcs.interpretAxis('time', var)
            print(var.axislist[idx].info())
            print('Time axis values:')
            print(timeax)
            print()
        except:
            print('Your variable probably does not have proper time axis.')

        # try get latitude axis
        try:
            latax=var.getLatitude()
            print('#'*50)
            print('Variable latitude axis:')
            print('#'*50)
            idx=funcs.interpretAxis('latitude', var)
            print(var.axislist[idx].info())
            print('Latitude axis values:')
            print(latax)
            print()
        except:
            print('Your variable probably does not have proper latitude axis.')

        # try get longitude axis
        try:
            lonax=var.getLongitude()
            print('#'*50)
            print('Variable longitude axis:')
            print('#'*50)
            idx=funcs.interpretAxis('longitude', var)
            print(var.axislist[idx].info())
            print('Longitude axis values:')
            print(lonax)
            print()
        except:
            print('Your variable probably does not have proper longitude axis.')

        try:
            units=var.units
        except:
            print('It is advised to write a "units" attribute to the data.')
        else:
            print('Data have unit of "%s"' %units)





