import cartopy.crs as ccrs

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
