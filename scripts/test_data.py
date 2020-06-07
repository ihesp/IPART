'''Check meta data of netcdf file.

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2020-03-31 10:42:41.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os, sys
import cdms2 as cdms
import numpy as np




#-------------Main---------------------------------
if __name__=='__main__':

    if len(sys.argv)<3:
        raise Exception("Please provide the path to the netcdf file, and the id of the variable.")

    file_path=os.path.abspath(sys.argv[1])

    if not os.path.exists(file_path):
        raise Exception("Input file not exists.")

    print('\n### <test_data>: Read in file:\n',file_path)

    try:
        fin=cdms.open(file_path,'r')
    except Exception as e:
        raise Exception("Failed to open file.\nException: %s" %str(e))

    try:
        var=fin[sys.argv[2]]
    except:
        raise Exception("Variable not found in file. Please verify the variable id is correct.")
    else:

        # check rank
        try:
            ndim=np.ndim(var)
            if ndim not in [3,4]:
                raise Exception("Data should be in rank 3 or 4, i.e. (time, lat, lon) or (time, level, lat, lon).")
        except:
            raise Exception("Data should be in rank 3 or 4, i.e. (time, lat, lon) or (time, level, lat, lon).")


        # try get time axis
        try:
            timeax=var.getTime()
            print('#'*50)
            print('Variable time axis:')
            print('#'*50)
            print(timeax.asComponentTime())
            print()
        except:
            print('Your variable probably does not have proper time axis.')

        # try get latitude axis
        try:
            latax=var.getLatitude()
            print('#'*50)
            print('Variable latitude axis:')
            print('#'*50)
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

    finally:
        fin.close()
    


    

