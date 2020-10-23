import numpy as np
import xarray as xr
import h5py

class GSFCmascons:
    def __init__(self, f, lon_wrap='pm180'):
        self.lat_centers = f['/mascon/lat_center'][0][:]
        self.lat_spans = f['/mascon/lat_span'][0][:]
        self.lon_centers = f['/mascon/lon_center'][0][:]
        self.lon_spans = f['/mascon/lon_span'][0][:]
        self.locations = f['/mascon/location'][0][:]
        self.basins = f['/mascon/basin'][0][:]
        self.areas = f['/mascon/area_km2'][0][:]
        self.cmwe = f['/solution/cmwe'][:]
        
        self.days_start = f['/time/ref_days_first'][0][:]
        self.days_middle = f['/time/ref_days_middle'][0][:]
        self.days_end = f['/time/ref_days_last'][0][:]
        self.times_start = self._set_times_as_datetimes(self.days_start)
        self.times_middle = self._set_times_as_datetimes(self.days_middle)
        self.times_end = self._set_times_as_datetimes(self.days_end)

        self.N_mascons = len(self.lat_centers)
        self.N_times = len(self.days_middle)
        self.labels = np.array([i for i in range(self.N_mascons)])
        
        self.reset_lon_bounds(lon_wrap)
        
        self.min_lats = self.lat_centers - self.lat_spans/2
        self.max_lats = self.lat_centers + self.lat_spans/2
        self.max_lats[self.min_lats < -90.0] = -89.5
        self.min_lats[self.min_lats < -90.0] = -90.0
        self.min_lats[self.max_lats > 90.0] = 89.5
        self.max_lats[self.max_lats > 90.0] = 90.0
        self.min_lons = self.lon_centers - self.lon_spans/2
        self.max_lons = self.lon_centers + self.lon_spans/2

    def reset_lon_bounds(self, lon_wrap):
        if lon_wrap == 'pm180':
            self.lon_centers[self.lon_centers > 180] -= 360
        elif lon_wrap == '0to360':
            self.lon_centers[self.lon_centers < 0] += 360

    def _set_times_as_datetimes(self, days):
        return np.datetime64('2002-01-01T00:00:00') + np.array([int(d*24) for d in days], dtype='timedelta64[h]')
    
    def as_dataset(self):
        ds = xr.Dataset({'cmwe': (['label', 'time'], self.cmwe),
                         'lat_centers': ('label', self.lat_centers),
                         'lat_spans': ('label', self.lat_spans),
                         'lon_centers': ('label', self.lon_centers),
                         'lat_spans': ('label', self.lon_spans),
                         'areas': ('label', self.areas),
                         'basins': ('label', self.basins),
                         'locations': ('label', self.locations),
                         'basins': ('label', self.basins),
                         'lats_max': ('label', self.max_lats),
                         'lats_min': ('label', self.min_lats),
                         'lons_max': ('label', self.max_lons),
                         'lons_min': ('label', self.min_lons),
                         'times_start': ('time', self.times_start),
                         'times_end': ('time', self.times_end),
                         'days_start': ('time', self.days_start),
                         'days_middle': ('time', self.days_middle),
                         'days_end': ('time', self.days_end)
                        }, coords={'label': self.labels, 'time': self.times_middle})
        return ds

def load_gsfc_solution(h5_filename, lon_wrap='pm180'):
    with h5py.File(h5_filename, mode='r') as f:
        mascons = GSFCmascons(f, lon_wrap)
    return mascons

def points_to_mascons(mascons, lats, lons, values):
    d2r = np.pi/180
    
    min_lats = mascons.lat_centers - mascons.lat_spans/2
    max_lats = mascons.lat_centers + mascons.lat_spans/2
    min_lons = mascons.lon_centers - mascons.lon_spans/2
    max_lons = mascons.lon_centers + mascons.lon_spans/2
    
    mscn_mean = np.nan * np.ones(mascons.N_mascons)
    for i in range(mascons.N_mascons):
        
        if np.min(lats) > max_lats[i]:
            continue
        if np.max(lats) < min_lats[i]:
            continue
        if np.min(lons) > max_lons[i]:
            continue
        if np.max(lons) < min_lons[i]:
            continue
        
        I_ = (lats >= min_lats[i]) & (lats < max_lats[i]) & (lons >= min_lons[i]) & (lons < max_lons[i])
        m = values[I_]
        m_lats = lats[I_]
        
        m_lats = m_lats[~np.isnan(m)]
        m = m[~np.isnan(m)]
        
        if len(m) == 0:
            continue
        if np.sum(~np.isnan(m)) == 0:
            continue

        cos_weight = np.cos(m_lats*d2r)
        mscn_mean[i] = np.nanmean(m) # * cos_weight) / (np.sum(cos_weight) * len(m))
    
    return mscn_mean

def calc_mascon_delta_cmwe(mascons, start_date, end_date):
    t_0 = np.datetime64(start_date)
    t_1 = np.datetime64(end_date)
    
    i_0 = np.abs(mascons.times_start - t_0).argmin()
    i_1 = np.abs(mascons.times_end - t_1).argmin()
    
    return mascons.cmwe[:,i_1] - mascons.cmwe[:,i_0]