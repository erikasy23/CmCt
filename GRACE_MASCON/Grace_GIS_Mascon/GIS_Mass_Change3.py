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

# Load GIS model into an Xarray
nc_filename = './lithk_GIS_ARC_PISM_asmb.nc'
gis_ds = xr.open_dataset(nc_filename, autoclose=True, engine='netcdf4')

lithk = gis_ds['lithk']
lithk_proj = gis_ds['Polar_Stereographic']

polar_stereographic = ccrs.Stereographic(
    central_latitude=90.0,
    central_longitude=-45.0,
    false_easting=0.0,
    false_northing=0.0,
    true_scale_latitude=70.0,
    globe=ccrs.Globe('WGS84')
)

# Compute mascon means
start_date = '2004-01-01' # 'YYYY-MM-DD'
end_date = '2016-01-01' # 'YYYY-MM-DD'
cmwe_delta = mascons.calc_mascon_delta_cmwe(gsfc, start_date, end_date)

# Select only GIS mascons
I_ = gsfc.locations == 1
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
ax.set_extent([-65, -20, 57, 84]) # Map bounds, [west, east, south, north]

# In order to set a colormap for the mascon plot, we first plot a scatter of the mascons.
# This will allow us to set colormap values for use with a second layer of "fill" data
# representing the total extent of each mascon. These scatter points will be covered by
# the "fill" images, and thus not visible in the final map.
sc = plt.scatter(lon_centers, lat_centers, 1, c=cmwe_delta, zorder=0, transform=ccrs.PlateCarree(),
                 cmap=plt.cm.RdBu, vmin=diverging_min, vmax=diverging_max)

normal = plt.Normalize(diverging_min, diverging_max)
cmap = plt.cm.RdBu(normal(cmwe_delta))

# Using the colormap info from above, draw each GIS mascon and fill with appropriate color
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
plt.title("GRACE Mascon Total Mass Change")
plt.savefig("GRACE_Mascon.png")

#plt.show()

# Transform projection to lat/lon
geodetic = ccrs.Geodetic(globe=ccrs.Globe('WGS84'))

yv, xv = np.meshgrid(gis_ds.y.data, gis_ds.x.data)

ll = geodetic.transform_points(src_crs=polar_stereographic, x=xv.flatten(), y=yv.flatten())
lons = ll[:,0]
lats = ll[:,1]

# # Calc mean field so we have a single field to plot
lithk_mean = lithk.mean(dim='time').data
lithk_mean = lithk_mean.transpose()
lithk_mean = lithk_mean.flatten()
#
# Calc difference between 2015.0 and 2000.0:
#lithk_delta = (lithk[3] - lithk[0]).data.transpose().flatten()

# Calc difference between 2015.0 and 2000.0:
lithk_delta = (lithk[3] - lithk[0]).data.transpose().flatten()

#start_date = '2004-01-01' # 'YYYY-MM-DD'
#end_date = '2014-01-01' # 'YYYY-MM-DD'
#
#lithk_start = lithk.interp(time=start_date).data.transpose().flatten()
#lithk_end = lithk.interp(time=end_date).data.transpose().flatten()
#
#lithk_delta = lithk_end - lithk_start

# Plot transformed data
fig = plt.figure(figsize=(5,6), dpi=300)

ax = plt.axes(projection=polar_stereographic)
ax.set_extent([-65, -20, 57, 84]) # Map bounds, [west, east, south, north]

sc = plt.scatter(lons, lats, 0.25, c=lithk_delta, transform=ccrs.Geodetic(), zorder=0, cmap=plt.cm.RdBu,
                 vmin=-8, vmax=8)

c = plt.colorbar(orientation='horizontal')
#c.set_label('Δ m ice, model ({0} to {1})'.format(start_date, end_date))

c.set_label('Model Delta thickness (m), 2004.0-2016.0')


ax.coastlines(resolution='50m', zorder=7, linewidth=0.5)
ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)
plt.title("Ice Sheet Model Total Mass Change")
plt.savefig("GIS_Model.png")

# Mascon-average lithk from GIS
lithk_delta[np.isnan(lithk_delta)] = 0
lithk_mascons = mascons.points_to_mascons(gsfc, lats, lons, lithk_delta)

# Ice thickness (m) to cm water equivalent:
rho_ice = 934 # kg/m^3
rho_water = 1000 # kg/m^3
lithk_mascons_cmwe = lithk_mascons * rho_ice / rho_water * 100

# I_ = ~np.isnan(lithk_mascons_cmwe)
I_ = gsfc.locations == 1
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

# Plot Mascon-Averaged GIS
plt.figure(figsize=(5,6), dpi=300)

ax = plt.axes(projection=polar_stereographic)
ax.set_extent([-65, -20, 57, 84]) # Map bounds, [west, east, south, north]

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
c.set_label('Model Delta cm water equiv., 2004.0-2016.0')

ax.coastlines(resolution='10m', zorder=7, linewidth=0.5)
ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)
plt.title("Ice Sheet Model Mascon Total Mass Change")
sc.remove()

plt.savefig("GIS_Model_mascon.png")
