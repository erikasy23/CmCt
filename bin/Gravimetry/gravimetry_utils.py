import os,sys
import glob as glob
import numpy as np
import h5py
import xarray as xr
import glob as glob
import requests
import re
from requests.exceptions import HTTPError, RequestException

import datetime

import mascons


import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')


# Method : Set model projection from standard definition
# Function to set projection based on selected loc
def set_projection(loc):
    if loc == "GIS":
        polar_stereographic = ccrs.Stereographic(
            central_latitude=90.0,
            central_longitude=-45.0,
            false_easting=0.0,
            false_northing=0.0,
            true_scale_latitude=70.0,
            globe=ccrs.Globe('WGS84')
        )
    elif loc == "AIS":
        polar_stereographic = ccrs.Stereographic(
            central_latitude=-90.0,
            central_longitude=0.0,
            false_easting=0.0,
            false_northing=0.0,
            true_scale_latitude=-71.0,
            globe=ccrs.Globe('WGS84')
        )
    return polar_stereographic



def loadGsfcMascons(Mascon_data_path):
           
    # Load GSFC mascons
    try:
        gsfc = mascons.load_gsfc_solution( Mascon_data_path, lon_wrap='pm180')
    except Exception as error:
        print('Error: Failed to load GSFC mascons.')
        print(error)
    return gsfc


def loadGisModel(nc_filename):
    # Load GIS model into an Xarray
    try:
        gis_ds = xr.open_dataset(nc_filename, autoclose=True, engine='netcdf4')
    except:
        print('Error: Failed to open model data; unexpected format found. Terminating calculation.')

    # check for the lithk variable
    try:
        lithk = gis_ds['lithk']
    except Exception as error:
        print('Error: lithk variable expected but not found in model. Terminating calculation.')
        print(error)

    return gis_ds


# Compute mascon means
def computeMasconMeans(gsfc, start_date, end_date, loc):

    global lat_centers
    global lon_centers
    global max_lons
    global min_lons
    global min_lats
    global max_lats
    global min_mscns
    global max_mscns
    global diverging_max
    global diverging_min
    global I_

    try:
        cmwe_delta = mascons.calc_mascon_delta_cmwe(gsfc, start_date, end_date)
    except Exception as error:
        print('Error: Failed to calculate mascon delta. Terminating calculation.')
        print(error)

    # Select only desired mascons
    if loc == "GIS":
        I_ = gsfc.locations == 1
    elif loc == "AIS":
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

    return cmwe_delta,I_



def transformToGeodetic(gsfc, gis_ds, start_date, end_date, polar_stereographic):
    # Put model into mascon space:

    # To compare with GRACE mascons, we need to compute lat/lon coordinates
    # for the grid locations and average them into the GSFC mascon boundaries.

    # First, we must transform from the original polar stereographic projection
    # into a geodetic lat/lon coordinate system. We plot the result of
    # this transformation to verify that the transformation was successful.

    # Then, we spatially average the data into mascon space and once more plot our result.

    # TODO: evaluate whether this transform has failed and return appropriate error

    # fetch the lithk variable from the model data structure
    lithk = gis_ds['lithk']

    # Transform projection to lat/lon
    geodetic = ccrs.Geodetic(globe=ccrs.Globe('WGS84'))

    yv, xv = np.meshgrid(gis_ds.y.data, gis_ds.x.data)

    ll = geodetic.transform_points(src_crs=polar_stereographic, x=xv.flatten(), y=yv.flatten())
    lons = ll[:,0]
    lats = ll[:,1]

    # # Calc difference between end_date and start_date:
    lithk_start = lithk.interp(time=start_date).data.transpose().flatten()
    lithk_end = lithk.interp(time=end_date).data.transpose().flatten()

    lithk_delta = lithk_end - lithk_start

    # Mascon-average lithk from GIS
    lithk_delta[np.isnan(lithk_delta)] = 0
    lithk_mascons = mascons.points_to_mascons(gsfc, lats, lons, lithk_delta)

    # Ice thickness (m) to cm water equivalent:
    rho_ice = 918 # kg/m^3
    rho_water = 1000 # kg/m^3
    lithk_mascons_cmwe = lithk_mascons * rho_ice / rho_water * 100

    # these variables depend only on the mascons here, which are fixed.
    mscns_trim = lithk_mascons_cmwe[I_]

    return mscns_trim, lithk_mascons_cmwe



def plotFigure(cmwe_delta, mscns_trim, cmwe_diff, gsfc, I_, lithk_mascons_cmwe,start_date, end_date, polar_stereographic,loc,plot_filename):

    
    current_directory = os.getcwd()
    parent_directory = os.path.abspath(os.path.join(current_directory, os.pardir, os.pardir))
    shapefile=os.path.join(parent_directory, 'data/ne_10m_coastline/ne_10m_coastline.shp')

    plt.figure(figsize=(24,14)) #, dpi=300)
    
    if loc == "GIS":
        extent = [-65, -20, 57, 84]
    elif loc == "AIS":
        extent = [-135, 45, -52, -52]
    else:
        print(f"Error: Input loc is equal to '{loc}', not 'AIS' or 'GIS")
        return None

    # Observed
    ax1 = plt.subplot(131, projection=polar_stereographic)
    ax1.set_extent(extent) # Map bounds, [west, east, south, north]

    sc = ax1.scatter(lon_centers, lat_centers, 1, c=cmwe_delta, zorder=0, transform=ccrs.PlateCarree(),
                     cmap=plt.cm.RdBu, vmin=diverging_min, vmax=diverging_max)

    normal = plt.Normalize(diverging_min, diverging_max)
    cmap = plt.cm.RdBu(normal(cmwe_delta))

    N_ints = 10
    for i in range(len(mscns_trim)):
        x = np.append(np.linspace(min_lons[i], max_lons[i], N_ints), np.linspace(max_lons[i], min_lons[i], N_ints))
        y = np.append(min_lats[i]*np.ones(N_ints), max_lats[i]*np.ones(N_ints))
        ax1.fill(x, y, facecolor=cmap[i][:], edgecolor='none', zorder=5, transform=ccrs.PlateCarree())

    c = plt.colorbar(sc, orientation='horizontal', ax=ax1, pad=0.04) #, fraction=0.046)
    c.set_label('cm water eq.', size=14)
    c.ax.tick_params(labelsize=12)

    ax1.add_geometries(list(shpreader.Reader(os.path.expanduser(shapefile)).geometries()), \
       ccrs.PlateCarree(), edgecolor='black', facecolor='none')
    # download coastline here: https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/

    # Add coastlines on top
    ax1.coastlines(resolution='10m', zorder=10, linewidth=0.5)
    # Add gridlines
    ax1.gridlines(zorder=5, linestyle=':', linewidth=0.5)

    ax1.set_title('Observed mass change\n({0} to {1})'.format(start_date, end_date), size=14)

    sc.remove()


    # Modeled
    ax2 = plt.subplot(132, projection=polar_stereographic)
    ax2.set_extent(extent) # Map bounds, [west, east, south, north]

    sc = ax2.scatter(gsfc.lon_centers[I_], gsfc.lat_centers[I_], 1, c=lithk_mascons_cmwe[I_], zorder=0,
                     transform=ccrs.PlateCarree(), cmap=plt.cm.RdBu, vmin=diverging_min, vmax=diverging_max)

    normal = plt.Normalize(diverging_min, diverging_max)
    cmap = plt.cm.RdBu(normal(lithk_mascons_cmwe[I_]))

    N_ints = 10
    for i in range(len(mscns_trim)):
        x = np.append(np.linspace(min_lons[i], max_lons[i], N_ints), np.linspace(max_lons[i], min_lons[i], N_ints))
        y = np.append(min_lats[i]*np.ones(N_ints), max_lats[i]*np.ones(N_ints))
        ax2.fill(x, y, facecolor=cmap[i][:], edgecolor='none', zorder=5, transform=ccrs.PlateCarree())

    c = plt.colorbar(sc, orientation='horizontal', ax=ax2, pad=0.04) #, fraction=0.046)
    c.set_label('cm water eq.', size=14)
    c.ax.tick_params(labelsize=12)

    ax2.add_geometries(list(shpreader.Reader(os.path.expanduser(shapefile)).geometries()), \
       ccrs.PlateCarree(), edgecolor='black', facecolor='none')

    # Add coastlines on top
    ax2.coastlines(resolution='10m', zorder=10, linewidth=0.5)
    # Add gridlines
    ax2.gridlines(zorder=5, linestyle=':', linewidth=0.5)

    # add model filename to subplot's title
    #ax2.set_title('Modeled mass change\n({0} to {1}\n{3})'.format(start_date, end_date, file_name[1]), size=14)
    ax2.set_title('Modeled mass change\n({0} to {1})'.format(start_date, end_date), size=14)

    sc.remove()

    # add some explanatory information
    plt.suptitle('Gravimetry Comparison Plots', fontsize=25)
    plt.subplots_adjust(top=0.83)


    # Obeserved-Modeled
    ax3 = plt.subplot(133, projection=polar_stereographic)
    ax3.set_extent(extent) # Map bounds, [west, east, south, north]

    sc = ax3.scatter(lon_centers, lat_centers, 1, c=cmwe_diff, zorder=0, transform=ccrs.PlateCarree(),
                     cmap=plt.cm.RdBu, vmin=diverging_min, vmax=diverging_max)

    normal = plt.Normalize(diverging_min, diverging_max)
    cmap = plt.cm.RdBu(normal(cmwe_diff))

    N_ints = 10
    for i in range(len(mscns_trim)):
        x = np.append(np.linspace(min_lons[i], max_lons[i], N_ints), np.linspace(max_lons[i], min_lons[i], N_ints))
        y = np.append(min_lats[i]*np.ones(N_ints), max_lats[i]*np.ones(N_ints))
        ax3.fill(x, y, facecolor=cmap[i][:], edgecolor='none', zorder=5, transform=ccrs.PlateCarree())

    c = plt.colorbar(sc, orientation='horizontal', ax=ax3, pad=0.04) #, fraction=0.046)
    c.set_label('cm water eq.', size=14)
    c.ax.tick_params(labelsize=12)

    ax3.add_geometries(list(shpreader.Reader(os.path.expanduser(shapefile)).geometries()), \
       ccrs.PlateCarree(), edgecolor='black', facecolor='none')

    # Add coastlines on top
    ax3.coastlines(resolution='10m', zorder=10, linewidth=0.5)
    # Add gridlines
    ax3.gridlines(zorder=5, linestyle=':', linewidth=0.5)

    # add model filename to subplot's title
    ax3.set_title('Residual mass change\n({0} to {1})'.format( start_date, end_date), size=14)

    sc.remove()

    # add some explanatory information
    plt.suptitle('Gravimetry Comparison Plots', fontsize=25)
    plt.subplots_adjust(top=1.93)#0.83
    
    plt.savefig(plot_filename)
    plt.show()

