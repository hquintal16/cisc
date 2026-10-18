#
# This script is for downloading ERA5 data from Copernicus
#
#   Author: Lauren Grimley, lauren.grimley@unc.edu
#   Last edited by: LEG 2/6/24
#                   HCQ 2024-05-06
#                   HCQ 2024-06-11
#
# Inputs:
#   1) required - the range of years to download the data for
#   2) required - bounding box
#   3) optional - variable names of interest (default is wind u-v, precip, temperature)
#   4) required - output directory
#
# Outputs:
#   netcdf files of era5 hourly data
#
# Dependencies:
#   Install the CDS API: https://cds.climate.copernicus.eu/api-how-to

import cdsapi
import os

# Edit: years, bounding box, working directory
year = list(range(1940, 2024)) #to get through 2023, write 2024 as last value
month = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
# bbox = [36.55, -84.35, 32.05, -75.45, ] 
bbox = [37, -85, 31, -75, ]#Carolinas
# bbox = [49.3457868,-125.000000,24.396308,-66.93457, ]#CONUS
out_dir = r'C:\wgetdown'

''' Run the script '''
# Loop through the years and download era5 using cdsapi.Client
for i in range(len(year)):
    for j in range(len(month)):
        file = str(year[i]) + '_' + str(month[j]) +'_era5.nc' #Carolinas
        fileout = os.path.join(out_dir, file)
        c = cdsapi.Client()
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'format': 'netcdf',
                'product_type':'reanalysis',
                # Update variables that you want to pull
                'variable': [
                    'total_precipitation',
                    '2m_dewpoint_temperature',
                    '2m_temperature',
                ],
                'year': [str(year[i]),
                        ],
                'month': [str(month[j]),
                         ],
                'day': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',
                    '25', '26', '27', '28', '29', '30', '31',
                ],
                'time': ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                         '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                         '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                         '18:00', '19:00', '20:00', '21:00', '22:00', '23:00',
                        ],
                'area': bbox,
            },
            fileout)
        print('Completed writing', file)
