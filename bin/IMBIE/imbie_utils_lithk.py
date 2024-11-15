import os
import numpy as np
import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

import cftime 
import datetime
from datetime import timedelta 




### Check the selected dates are within the range of model data
def check_datarange(time_var,start_date, end_date):

    calendar_type = time_var.to_index().calendar
    
    start_date_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_date_dt = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    
    # Adjust day to be 30 ( to avoid error if it's the 31st day in a 360_day calendar)
    start_date_cftime = cftime.datetime(start_date_dt.year, start_date_dt.month, min(start_date_dt.day, 30), calendar=calendar_type)    
    
    end_date_cftime = cftime.datetime(end_date_dt.year, end_date_dt.month, min(end_date_dt.day, 30), calendar=calendar_type)    

    
    # Get the minimum and maximum values directly from the time variable
    min_time = time_var.values.min()
    max_time = time_var.values.max()
    
    # Check if the selected start and end dates are within the range
    if min_time <= start_date_cftime <= max_time and min_time <= end_date_cftime <= max_time:
        print(f"The selected dates {start_date} and {end_date} are within the range of the model data.")
    else:
        raise ValueError(f"Error: The selected dates {start_date} or {end_date} are out of range. Model data time range is from {min_time} to {max_time}.")
    


### Load the model data and calculate  model mass balance for each basin and total mass balance for whole region
## Interpolate the data to each imbie time and calculate the time varying mass change
def process_model_data(nc_filename,IMBIE_total_mass_change_sum,start_date, end_date,rho_ice,projection,shape_filename,icesheet):
    #Model data
    gis_ds = xr.open_dataset(nc_filename)
    lithk = gis_ds['lithk']
    time_var = gis_ds['time']

    # Load basin shapefile 
    basins_gdf = gpd.read_file(shape_filename)
    
    # Check the selcted dates are within the range of model data
    check_datarange(time_var,start_date, end_date)
        
     # Set start_date as the first date in 'Year' and filtered_time_var as all subsequent dates
    start_date_imbie = IMBIE_total_mass_change_sum['Year'].iloc[0]
    filtered_time_var = IMBIE_total_mass_change_sum['Year'].iloc[1:]
    
    #Chnage to modal data format
    start_date_dt_imbie = datetime.datetime.strptime(start_date_imbie, '%Y-%m-%d')   
    # Adjust day to be 30 ( to avoid error if it's the 31st day in a 360_day calendar)
    start_date_cftime = cftime.datetime(start_date_dt_imbie.year, start_date_dt_imbie.month, min(start_date_dt_imbie.day, 30), calendar=gis_ds.time.to_index().calendar)    
    
    
    
    #calculate area = x_resolution*y_resolution
    x_coords = gis_ds['x'].values
    y_coords = gis_ds['y'].values
    x_resolution = abs(x_coords[1] - x_coords[0])
    y_resolution = abs(y_coords[1] - y_coords[0])
    
    # Create a list of Point geometries from coordinate grids
    points = [Point(x, y) for x in x_coords for y in y_coords]
    
    # Initialize a dictionary to store residuals
    modal_mass_changes = {}
    
    # Interpolate limnsw at the start date (initial reference)
    lithk_start = lithk.interp(time=start_date_cftime).data.transpose().flatten()
    
    
    # Loop through each filtered time step to calculate the residual
    for i, time_step in enumerate(filtered_time_var):
    
        # Parse the current time step as a datetime object
        time_step_imbie = datetime.datetime.strptime(time_step, '%Y-%m-%d')
    
        # Convert the parsed time step to a cftime.datetime with the correct calendar
        # Adjust day to be 30 ( to avoid error if it's the 31st day in a 360_day calendar)
        time_step_cftime = cftime.datetime(time_step_imbie.year, time_step_imbie.month, min(time_step_imbie.day, 30), 
                                           calendar=gis_ds.time.to_index().calendar)  
     
        # Interpolate limnsw at the current time_step
        lithk_current = lithk.interp(time=time_step_cftime).data.transpose().flatten()
        
        # Calculate the residual (difference from start)
        lithk_delta = lithk_current - lithk_start
        lithk_start=lithk_current
    
        lithk_delta[np.isnan(lithk_delta)] = 0
        
        #calculate area = x_resolution*y_resolution
        lithk_delta = (lithk_delta * x_resolution*y_resolution)*rho_ice * 1e-12
        
        lithk_delta_flat = lithk_delta.flatten()
        
        # Create a lithk_df DataFrame with Geometry and Values
        lithk_df = pd.DataFrame({
            'geometry': points,
            'lithk_delta': lithk_delta_flat
        })
    
        # Convert lithk_df DataFrame to lithk_gdf GeoDataFrame
        lithk_gdf = gpd.GeoDataFrame(lithk_df, geometry='geometry', crs=projection)    
        
    
           # Perform the spatial join only once in the first iteration
        if i == 0:
            joined_gdf = gpd.sjoin(lithk_gdf, basins_gdf, how="inner", predicate='intersects')
    
    
        # Update the lithk_delta in joined_gdf
        joined_gdf['lithk_delta'] = lithk_gdf['lithk_delta']
           
        # Sum lithk_delta values by basin
        if icesheet == "GIS":
             # Sum lithk_delta values by subregion column
            basin_mass_change_sums = joined_gdf.groupby('SUBREGION1')['lithk_delta'].sum()
            # Sum lithk_delta values by the 'Regions' column
            region_mass_change_sums = None  # No regions for Greenland
        elif icesheet == "AIS":
            # Sum lithk_delta values by subregion column
            basin_mass_change_sums = joined_gdf.groupby('Subregion')['lithk_delta'].sum()
            # Sum lithk_delta values by the 'Regions' column
            region_mass_change_sums = joined_gdf.groupby('Regions')['lithk_delta'].sum()
        else:
            raise ValueError("Invalid iceshee value. Must be 'GIS' or 'AIS'.")
        
        # Sum all of the basin mass change
        model_total_mass_balance= basin_mass_change_sums.sum()
          
        # Store the residual in the dictionary with the time as the key
        modal_mass_changes[str(time_step)] = model_total_mass_balance
        
    # Return results
    return modal_mass_changes





### IMBIE data date format conversion
# Define a function to convert fractional years to a precise datetime format
def fractional_year_to_date(year):
    year_int = int(year)  # Extract the integer part (the full year)
    fraction = year - year_int  # Extract the fractional part
    
    # Start at the beginning of the year
    start_of_year = pd.Timestamp(f'{year_int}-01-01')
    
    # Determine if it's a leap year
    if pd.Timestamp(f'{year_int}-12-31').is_leap_year:
        total_days_in_year = 366
    else:
        total_days_in_year = 365
    
    # Convert the fractional part into the corresponding number of days
    fractional_days = fraction * total_days_in_year
    
    # Add the fractional days to the start of the year to get the correct date
    return start_of_year + timedelta(days=fractional_days)


# Group the data by year
def assign_month_order(group):
    # Get the month of the first entry for the year
    first_month = group['Date'].dt.month.iloc[0]
    
    # Create a month order starting from the first month and increasing by 1 for each subsequent entry
    group['Month_Order'] = range(first_month, first_month + len(group))
    return group



### Extract time varying IMBIE mass balance data and calculate the time varying mass difference 
def sum_MassBalance(obs_filename,start_date,end_date,mass_balance_column):
    
    # Load the CSV file
    mass_balance_data = pd.read_csv(obs_filename)
    
    # Column names
    date_column = 'Year'
    
    # Ensure the 'Year' column is treated as float to capture the fractional year part
    mass_balance_data['Year'] = mass_balance_data['Year'].astype(float)
    
    # Apply the conversion function to the 'Year' column
    mass_balance_data['Date'] = mass_balance_data['Year'].apply(fractional_year_to_date)
    
    # Sort the data by 'Date' column to ensure it’s in increasing order of both year and fraction
    mass_balance_data = mass_balance_data.sort_values(by='Date')
      
    # Apply the function to each group of data (grouped by the year)
    mass_balance_data = mass_balance_data.groupby(mass_balance_data['Date'].dt.year).apply(assign_month_order)

    
    # Convert 'Year' column to year-month-day format where month is 'Month_Order'    
    mass_balance_data['Year'] = mass_balance_data.apply(lambda row: f"{row['Date'].year}-{str(row['Date'].month).zfill(2)}-{str(row['Date'].day).zfill(2)}", axis=1) 
    # mass_balance_data['Year'] = mass_balance_data.apply(lambda row: f"{row['Date'].year}-{str(row['Month_Order']).zfill(2)}-01", axis=1)    
    
    # Reset the index to flatten the multi-index structure
    mass_balance_data = mass_balance_data.reset_index(drop=True)
    
    
    # Check if the column exists in the DataFrame
    if mass_balance_column not in mass_balance_data.columns:
        raise ValueError(f"Error: The column '{mass_balance_column}' does not exist in the CSV file.")
    
    
    # Get the initial mass balance value before the start date
    data_before_start_date = mass_balance_data[mass_balance_data['Year'] < start_date]
    if data_before_start_date.empty:
        raise ValueError(f"Error: No data available before the start date {start_date}.")
    mass_balance_start_value = data_before_start_date[mass_balance_column].iloc[-1]  # Last value before start date
    
    # Filter data between start_date_converted and end_date_converted (inclusive)
    filtered_data = mass_balance_data[
        (mass_balance_data['Year'] >= start_date) & (mass_balance_data['Year'] <= end_date)
    ]


    # Initialize the previous day's mass balance value to the starting mass balance
    previous_mass_balance = mass_balance_start_value
    
    # Calculate daily mass change from the previous day for each time step
    mass_changes = []  # To store the daily mass changes
    
    for index, row in filtered_data.iterrows():
        current_mass_balance = row[mass_balance_column]
        # Calculate the change from the previous day's balance
        mass_change = current_mass_balance - previous_mass_balance
        mass_changes.append(mass_change)
        # Update previous_mass_balance for the next iteration
        previous_mass_balance = current_mass_balance
    
    # Assign the calculated mass changes to a new column in the DataFrame
    filtered_data['IMBIE_Mass_Change'] = mass_changes

    return filtered_data







