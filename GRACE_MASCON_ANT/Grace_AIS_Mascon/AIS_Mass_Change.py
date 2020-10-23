#from ctypes import *
#lib1 = cdll.LoadLibrary('/opt/python3/lib')
#
import os

#os.environ['LD_LIBRARY_PATH'] = '/opt/python3/lib'

#import matplotlib as mpl
#if os.environ.get('DISPLAY','') == '':
#    print('no display found. Using non-interactive Agg backend')
#    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mascons

# Load GSFC mascons
h5_filename = './GSFC.glb.200301_201607_v02.4-GeruoA.h5'
gsfc = mascons.load_gsfc_solution(h5_filename, lon_wrap='pm180')

# Load AIS model into an Xarray
#   NOTE: Need to set decode_times=False because invalid times in NetCDF
#     "ValueError: unable to decode time units 'seconds since 1990-01-01 00:00:00'
#      with calendar 'standard'. Try opening your dataset with decode_times=False."
# nc_filename = './lithk_AIS_PSU_EQMEC_asmb.nc'#lithk_AIS_ARC_PISM1_asmb.nc'
nc_filename = './lithk_AIS_PSU_EQMEC_asmb.nc'#lithk_AIS_ARC_PISM1_asmb.nc'
ais_ds = xr.open_dataset(nc_filename, autoclose=True, engine='netcdf4', decode_times=False)

# NOTE: Include this workaround for invalid times in AIS files:
times = ['2000-01-01T00:00:00.000000000','2004-12-31T05:03:50.000000000','2009-12-31T10:07:40.000000000',
         '2014-12-31T15:11:30.000000000','2019-12-31T20:15:20.000000000','2024-12-31T01:19:10.000000000',
         '2029-12-31T06:23:00.000000000','2034-12-31T11:26:50.000000000','2039-12-31T16:30:40.000000000',
         '2044-12-30T21:34:30.000000000','2049-12-31T02:38:20.000000000','2054-12-31T07:42:10.000000000',
         '2059-12-31T12:46:00.000000000','2064-12-30T17:49:50.000000000','2069-12-30T22:53:40.000000000',
         '2074-12-31T03:57:30.000000000','2079-12-31T09:01:20.000000000','2084-12-30T14:05:10.000000000',
         '2089-12-30T19:09:00.000000000','2094-12-31T00:12:50.000000000','2099-12-31T05:16:40.000000000']
ais_ds['time'] = [np.datetime64(t) for t in times]


lithk = ais_ds['lithk']
#lithk_proj = ais_ds['mapping']

#Set Polar Sterographic Projection definition:
#Projection can be defined from the loaded model or by setting the definition independent of the model. Since the CmCt uses a standard projection, it is probably best to use the second method.

# # Method 1: Set model projection from model projection information
# polar_stereographic = ccrs.Stereographic(
#     central_latitude=lithk_proj.latitude_of_projection_origin,
#     central_longitude=lithk_proj.straight_vertical_longitude_from_pole,
#     false_easting=lithk_proj.false_easting,
#     false_northing=lithk_proj.false_northing,
#     true_scale_latitude=lithk_proj.standard_parallel,
#     globe=ccrs.Globe('WGS84')
# )

# Method 2: Set model projection from standard definition
polar_stereographic = ccrs.Stereographic(
                                         central_latitude=-90.0,
                                         central_longitude=0.0,
                                         false_easting=0.0,
                                         false_northing=0.0,
                                         true_scale_latitude=-71.0,
                                         globe=ccrs.Globe('WGS84')
                                         )

# Compute mascon means
start_date = '2004-01-01' # 'YYYY-MM-DD'
end_date = '2016-01-01' # 'YYYY-MM-DD'
cmwe_delta = mascons.calc_mascon_delta_cmwe(gsfc, start_date, end_date)

# Compute mascon means
start_date = '2004-01-01' # 'YYYY-MM-DD'
end_date = '2014-01-01' # 'YYYY-MM-DD'
cmwe_delta = mascons.calc_mascon_delta_cmwe(gsfc, start_date, end_date)

# Select only AIS mascons
I_ = (gsfc.locations == 3) | (gsfc.locations == 4)
cmwe_delta = cmwe_delta[I_]
lat_centers = gsfc.lat_centers[I_]
lon_centers = gsfc.lon_centers[I_]
min_lons = gsfc.min_lons[I_]
max_lons = gsfc.max_lons[I_]
min_lats = gsfc.min_lats[I_]
max_lats = gsfc.max_lats[I_]

min_mscns = np.min(cmwe_delta)
max_mscns = np.max(cmwe_delta)

diverging_max = np.max([np.abs(min_mscns), np.abs(max_mscns)])
diverging_min = -diverging_max

# Create figure and set projection:
plt.figure(figsize=(5,5), dpi=300)

ax = plt.axes(projection=polar_stereographic)

# In order to set a colormap for the mascon plot, we first plot a scatter of the mascons.
# This will allow us to set colormap values for use with a second layer of "fill" data
# representing the total extent of each mascon. These scatter points will be covered by
# the "fill" images, and thus not visible in the final map.
sc = plt.scatter(lon_centers, lat_centers, 1, c=cmwe_delta, zorder=0, transform=ccrs.PlateCarree(),
                 cmap=plt.cm.RdBu, vmin=diverging_min, vmax=diverging_max)

normal = plt.Normalize(diverging_min, diverging_max)
cmap = plt.cm.RdBu(normal(cmwe_delta))

# Using the colormap info from above, draw each AIS mascon and fill with appropriate color
N_ints = 10
for i in range(len(cmwe_delta)):
    x = np.append(np.linspace(min_lons[i], max_lons[i], N_ints),
                  np.linspace(max_lons[i], min_lons[i], N_ints))
    y = np.append(min_lats[i]*np.ones(N_ints), max_lats[i]*np.ones(N_ints))
    plt.fill(x, y, facecolor=cmap[i][:], edgecolor='none', zorder=5, transform=ccrs.PlateCarree())

c = plt.colorbar(orientation='horizontal')
#c.set_label('Δ cm water equiv., GRACE ({0} to {1})'.format(start_date, end_date))
c.set_label('GRACE Delta cm water equiv., 2004.0-2016.0')

ax.coastlines(resolution='50m', zorder=7, linewidth=0.5)
ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)
plt.title("GRACE Mascon Total Mass Change ANT")
plt.savefig("GRACE_Mascon_ANT.png")

#plt.show()

# Transform projection to lat/lon
geodetic = ccrs.Geodetic(globe=ccrs.Globe('WGS84'))

yv, xv = np.meshgrid(ais_ds.y.data, ais_ds.x.data)

ll = geodetic.transform_points(src_crs=polar_stereographic, x=xv.flatten(), y=yv.flatten())
lons = ll[:,0]
lats = ll[:,1]

start_date = '2004-01-01' # 'YYYY-MM-DD'
end_date = '2014-01-01' # 'YYYY-MM-DD'

lithk_start = lithk.interp(time=start_date).data.transpose().flatten()
lithk_end = lithk.interp(time=end_date).data.transpose().flatten()
lithk_delta = lithk_end - lithk_start

# Plot transformed data
fig = plt.figure(figsize=(5,6), dpi=300)

ax = plt.axes(projection=polar_stereographic)

sc = plt.scatter(lons, lats, 0.25, c=lithk_delta, transform=ccrs.Geodetic(), zorder=0, cmap=plt.cm.RdBu,
                 vmin=-4, vmax=4)

c = plt.colorbar(orientation='horizontal')
c.set_label('Δ m ice, model ({0} to {1})'.format(start_date, end_date))

ax.coastlines(resolution='50m', zorder=7, linewidth=0.5)
ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)
plt.title("Ice Sheet Model Total Mass Change ANT")
plt.savefig("AIS_Model.png")

# Mascon-average lithk from AIS
lithk_delta[np.isnan(lithk_delta)] = 0
lithk_mascons = mascons.points_to_mascons(gsfc, lats, lons, lithk_delta)

# Ice thickness (m) to cm water equivalent:
rho_ice = 934 # kg/m^3
rho_water = 1000 # kg/m^3
lithk_mascons_cmwe = lithk_mascons * rho_ice / rho_water * 100


#Plot Ice Sheet Model for Antarctica
I_ = (gsfc.locations == 3) | (gsfc.locations == 4)

mscns_trim = lithk_mascons_cmwe[I_]
min_lats = gsfc.min_lats[I_]
max_lats = gsfc.max_lats[I_]
min_lons = gsfc.min_lons[I_]
max_lons = gsfc.max_lons[I_]
# min_mscns = np.min(mscns_trim)
# max_mscns = np.max(mscns_trim)
min_mscns = np.min(cmwe_delta) # min/max from Mascons to have comparable plot.
max_mscns = np.max(cmwe_delta)

diverging_max = np.max([np.abs(min_mscns), np.abs(max_mscns)])
diverging_min = -diverging_max

# Plot Mascon-Averaged AIS
plt.figure(figsize=(5,6), dpi=300)

ax = plt.axes(projection=polar_stereographic)

sc = plt.scatter(gsfc.lon_centers, gsfc.lat_centers, 1, c=lithk_mascons_cmwe, zorder=0,
                 transform=ccrs.PlateCarree(), cmap=plt.cm.RdBu, vmin=diverging_min, vmax=diverging_max)

normal = plt.Normalize(diverging_min, diverging_max)
cmap = plt.cm.RdBu(normal(mscns_trim))

N_ints = 10
for i in range(len(mscns_trim)):
    x = np.append(np.linspace(min_lons[i], max_lons[i], N_ints), np.linspace(max_lons[i], min_lons[i], N_ints))
    y = np.append(min_lats[i]*np.ones(N_ints), max_lats[i]*np.ones(N_ints))
    plt.fill(x, y, facecolor=cmap[i][:], edgecolor='none', zorder=5, transform=ccrs.PlateCarree())

c = plt.colorbar(orientation='horizontal')
c.set_label('Δ cm w.e., model ({0} to {1})'.format(start_date, end_date))



ax.coastlines(resolution='10m', zorder=7, linewidth=0.5)
ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)
plt.title("Ice Sheet Model Mascon Total Mass Change ANT")
sc.remove()

plt.savefig("AIS_Model_mascon.png")
