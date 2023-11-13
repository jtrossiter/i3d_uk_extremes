"""
Gets the area averaged june daily temperatures for each day of each simulation in a given ensemble, then calculates the daily average, max,
and min across all simulations, and also the monthly averages for each simulation
"""

import numpy as np
import xarray as xr
import os

def get_temperature_data(file):
    dataset = xr.load_dataset(file)
    # load landfrac file
    lf = np.load('landfrac.npy')
    lf = xr.DataArray(lf)
    lf = lf.rename({'dim_0':'lat', 'dim_1':'lon'})
    weights = lf
    weights.name = 'weights'

    # UK land mass bounding box
    lat_min = 49.959999905
    lat_max = 58.6350001085
    lon_min = 360-7.57216793459
    lon_max = 1.68153079591

    # find nearest grid points to each of these lats/longs
    lats = dataset.lat.to_numpy()
    longs = dataset.lon.to_numpy()

    lat_min = lats[np.abs(lats - lat_min).argmin()]
    lat_max = lats[np.abs(lats - lat_max).argmin()]
    lon_min = longs[np.abs(longs - lon_min).argmin()]
    lon_max = longs[np.abs(longs - lon_max).argmin()]

    # slice along relevant time, lat, and lon
    dataset = dataset.sel(lat=slice(lat_min,lat_max))
    dataset = dataset.sel(time='2023-06') # only june 2023
    dataset = xr.concat([dataset.sel(lon=slice(lon_min,357.5)), dataset.sel(lon=slice(0,lon_max))], dim="lon")

    temps = np.zeros(30)
    for i in range(30):
         # area averaged temperature across each day weighted by landfrac
        dataset_weighted = dataset.TREFHTMX.sel(time = '2023-06-'+str(i+1).zfill(2)+' 00:00:00').weighted(weights)
        temps[i] = dataset_weighted.mean(dim=['lat','lon'])
    return list(temps), [np.mean(temps)]
"""
dirs = ['Historical/s1946411','Historical/jubauer', 'Historical/s1935349']
hist_daily_june_all = []
hist_monthly_avg = []
i = 0
for dir in dirs:
    for filename in os.listdir(dir):
        f = os.path.join(dir, filename)
        hist_daily_june_all = hist_daily_june_all + get_temperature_data(f)[0]
        hist_monthly_avg = hist_monthly_avg + get_temperature_data(f)[1]
        if (i == 0):
            hist_daily_june_avg = np.array(get_temperature_data(f)[0])
        else:
            hist_daily_june_avg += np.array(get_temperature_data(f)[0])
        i += 1
hist_daily_june_avg = hist_daily_june_avg/75

dirs = ['Natural/s1946411','Natural/eleanorsenior', 'Natural/krismakay12']
nat_daily_june_all = []
nat_monthly_avg = []
i = 0
for dir in dirs:
    for filename in os.listdir(dir):
        f = os.path.join(dir, filename)
        nat_daily_june_all = nat_daily_june_all + get_temperature_data(f)[0]
        nat_monthly_avg = nat_monthly_avg + get_temperature_data(f)[1]
        if (i == 0):
            nat_daily_june_avg = np.array(get_temperature_data(f)[0])
        else:
            nat_daily_june_avg += np.array(get_temperature_data(f)[0])
        i += 1
nat_daily_june_avg = nat_daily_june_avg/75

hist_daily_june_all = np.array(hist_daily_june_all) - 273.15 - 1 # convert to kelvin and account for bias
hist_monthly_avg = np.array(hist_monthly_avg) - 273.15 - 1
nat_daily_june_all = np.array(nat_daily_june_all) - 273.15 - 1
nat_monthly_avg = np.array(nat_monthly_avg) - 273.15 - 1

np.save('hist_daily_june_all', hist_daily_june_all)
np.save('nat_daily_june_all', nat_daily_june_all)
np.save('hist_monthly_avg', hist_monthly_avg)
np.save('nat_monthly_avg', nat_monthly_avg)
"""
# get daily average, max, and min temps
hist_daily_june_avg = np.zeros(30)
hist_daily_june_max = np.zeros(30)
hist_daily_june_min = np.zeros(30)
nat_daily_june_avg = np.zeros(30)
nat_daily_june_max = np.zeros(30)
nat_daily_june_min = np.zeros(30)

hist_daily_june_all = np.load('hist_daily_june_all.npy')
nat_daily_june_all = np.load('nat_daily_june_all.npy')
for i in range(30):
    temp1 = np.zeros(30)
    temp2 = np.zeros(30)
    for j in range(30):
        # get temp values for each day
        temp1[j] = hist_daily_june_all[j*30+i]
        temp2[j] = nat_daily_june_all[j*30+i]
    hist_daily_june_avg[i] = np.mean(temp1)
    hist_daily_june_max[i] = np.max(temp1)
    hist_daily_june_min[i] = np.min(temp1)
    nat_daily_june_avg[i] = np.mean(temp2)
    nat_daily_june_max[i] = np.max(temp2)
    nat_daily_june_min[i] = np.min(temp2)

np.save('hist_daily_june_avg', hist_daily_june_avg)
np.save('hist_daily_june_max', hist_daily_june_max)
np.save('hist_daily_june_min', hist_daily_june_min)

np.save('nat_daily_june_avg', nat_daily_june_avg)
np.save('nat_daily_june_max', nat_daily_june_max)
np.save('nat_daily_june_min', nat_daily_june_min)