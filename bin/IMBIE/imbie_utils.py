import os
import numpy as np
import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

import cftime 
import datetime
from datetime import timedelta 



### #Adjust the start and end date according to the time variable of model data
def adjust_for_calendar(calendar_type, date_dt, time_var):
    """
    Adjusts the date format to match the type of the time variable in the dataset.
    If the calendar is '360_day', adjust the day component accordingly.
    """
    # If the time_var is in numpy.datetime64, return the date as np.datetime64
    if isinstance(time_var.values[0], np.datetime64):
        return np.datetime64(date_dt)

    # If the time_var uses a cftime calendar, handle accordingly
    elif isinstance(time_var.values[0], cftime.datetime):
        if calendar_type == '360_day' and date_dt.day > 30:
            # In a 360-day calendar, each month has only 30 days.
            return cftime.datetime(date_dt.year, date_dt.month, 30, calendar=calendar_type)
        else:
            # For other calendars, use the original date.
            return cftime.datetime(date_dt.year, date_dt.month, date_dt.day, calendar=calendar_type)

    else:
        raise TypeError("Unsupported time format in the dataset.")

def convert_to_model_calendar(time_var, start_date, end_date):
    # Parse input dates to datetime objects
    start_date_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_date_dt = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    
    # Extract the calendar type used by the time variable
    try:
        calendar_type = time_var.to_index().calendar
        # print(calendar_type)
    except AttributeError:
        # Default to the 'standard' or 'gregorian' calendar if calendar attribute doesn't exist
        calendar_type = 'standard'
    
    # Convert the start_date and end_date to the correct calendar type
    start_date_cftime = adjust_for_calendar(calendar_type, start_date_dt,time_var)
    end_date_cftime = adjust_for_calendar(calendar_type, end_date_dt,time_var)
    
    return start_date_cftime, end_date_cftime

### Check the selected dates are within the range of model data
def check_datarange(time_var,start_date, end_date):

    # Get the minimum and maximum values directly from the time variable
    min_time = time_var.values.min()
    max_time = time_var.values.max()
    
    fomatted_start_date, fomatted_end_date =convert_to_model_calendar(time_var, start_date, end_date)
    
    # Check if the selected start and end dates are within the range
    if min_time <= fomatted_start_date <= max_time and min_time <= fomatted_end_date <= max_time:
        print(f"The selected dates {start_date} and {end_date} are within the range of the model data.")
    else:
        raise ValueError(f"Error: The selected dates {start_date} or {end_date} are out of range. Model data time range is from {min_time} to {max_time}.")




### Load the model data and calculate  model mass balance for each basin and total mass balance for whole region
def process_model_data(nc_filename,start_date, end_date,rho_ice,projection,shape_filename,icesheet):
    #Model data
    gis_ds = xr.open_dataset(nc_filename)
    lithk = gis_ds['lithk']
    time_var = gis_ds['time']
    
    # Check the selcted dates are within the range of model data
    check_datarange(time_var,start_date, end_date)
        
    # Interpolate lithk values at the start and end dates
    lithk_start = lithk.interp(time=start_date).data.transpose().flatten()
    lithk_end = lithk.interp(time=end_date).data.transpose().flatten()
    
    # Calculate the difference
    lithk_delta = lithk_end - lithk_start
    
    # Replace NaN values with 0
    lithk_delta[np.isnan(lithk_delta)] = 0
    
    
    # Change Ice thickness unit from (m) to mass (kg) to gigatonnes(Gt)
    # ice thickness*area* density of ice* 1e-12
    #calculate area = x_resolution*y_resolution
    x_coords = gis_ds['x'].values
    y_coords = gis_ds['y'].values
    x_resolution = abs(x_coords[1] - x_coords[0])
    y_resolution = abs(y_coords[1] - y_coords[0])
    
    lithk_delta = (lithk_delta * x_resolution*y_resolution)*rho_ice * 1e-12
    
    
    # Create a list of Point geometries from coordinate grids
    points = [Point(x, y) for x in x_coords for y in y_coords]
    
    # Flatten lithk_delta to match the points list 
    lithk_delta_flat = lithk_delta.flatten()
    
    # Create DataFrame
    lithk_df = pd.DataFrame({
        'geometry': points,
        'lithk_delta': lithk_delta_flat
    })
    
    # Convert DataFrame to GeoDataFrame
    lithk_gdf = gpd.GeoDataFrame(lithk_df, geometry='geometry', crs=projection)
    
    # Load basin shapefile 
    basins_gdf = gpd.read_file(shape_filename)
    
    # Perform spatial join
    joined_gdf = gpd.sjoin(lithk_gdf, basins_gdf, how="inner", predicate='intersects')
    
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
    
    # Return all results as a dictionary
    return {
        'model_total_mass_balance': model_total_mass_balance,
        'basin_mass_change_sums': basin_mass_change_sums,
        'region_mass_change_sums': region_mass_change_sums
    }






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



### Extract IMBIE mass balance data
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
    
    # Convert 'Year' column to year-month-01 format where month is 'Month_Order'
    mass_balance_data['Year'] = mass_balance_data.apply(lambda row: f"{row['Date'].year}-{str(row['Month_Order']).zfill(2)}-01", axis=1)
    
    # Reset the index to flatten the multi-index structure
    mass_balance_data = mass_balance_data.reset_index(drop=True)

    
    # Check if the column exists in the DataFrame
    if mass_balance_column not in mass_balance_data.columns:
        raise ValueError(f"Error: The column '{mass_balance_column}' does not exist in the CSV file.")

    
    # Filter the data for the end date
    end_data = mass_balance_data[mass_balance_data['Year'] == end_date]    
    if end_data.empty:
        raise ValueError(f"Error: No data available for the end date {end_date}.")
    mass_balance_end_value = end_data[mass_balance_column].iloc[-1]  # Last value before or at the end date

    
    # Filter the data for one date before the start date
    data_before_start_date = mass_balance_data[mass_balance_data[date_column] < start_date]
    if data_before_start_date.empty:
        raise ValueError(f"Error: No data available before the start date {start_date}.")
    mass_balance_start_value = data_before_start_date[mass_balance_column].iloc[-1]  # Last value before start date
    
    # Subtract the two values to get the total mass balance change
    IMBIE_total_mass_change_sum = mass_balance_end_value - mass_balance_start_value
    
    return IMBIE_total_mass_change_sum



### Calculate mass balance difference of IMBIE and model data
def process_IMBIE(obs_filename, start_date, end_date, icesheet, basin_result,mass_balance_column,obs_east_filename=None, obs_west_filename=None, obs_peninsula_filename=None):
    results = {}

    # model mass balance
    model_total_mass_balance = basin_result['model_total_mass_balance']
    
    # IMBIE total mass balance
    IMBIE_total_mass_change_sum = sum_MassBalance(obs_filename, start_date, end_date,mass_balance_column)
    
    # Calculate difference of IMBIE-model mass change
    delta_masschange = IMBIE_total_mass_change_sum - model_total_mass_balance
    
    # Store total mass balance results in the dictionary
    results['IMBIE_total_mass_change_sum'] = IMBIE_total_mass_change_sum
    results['delta_masschange'] = delta_masschange
    
    # Check if all required (regional) files are available for Antarctica
    if icesheet == "AIS":
        region_mass_change_sums = basin_result.get('region_mass_change_sums') 
        print_regionalresult_check = 'NO'
        if (obs_east_filename and os.path.exists(obs_east_filename)) and \
           (obs_west_filename and os.path.exists(obs_west_filename)) and \
           (obs_peninsula_filename and os.path.exists(obs_peninsula_filename)):
            
            print_regionalresult_check = 'YES' 
            
            # Calculate total mass for each region
            IMBIE_total_mass_change_sum_east = sum_MassBalance(obs_east_filename, start_date,         end_date,mass_balance_column)
            IMBIE_total_mass_change_sum_west = sum_MassBalance(obs_west_filename, start_date, end_date,mass_balance_column)
            IMBIE_total_mass_change_sum_peninsula = sum_MassBalance(obs_peninsula_filename, start_date, end_date,mass_balance_column)
            
            # Calculate the difference of IMBIE-model mass change for each region
            delta_masschange_east = IMBIE_total_mass_change_sum_east - region_mass_change_sums['East']
            delta_masschange_west = IMBIE_total_mass_change_sum_west - region_mass_change_sums['West']
            delta_masschange_peninsula = IMBIE_total_mass_change_sum_peninsula - region_mass_change_sums['Peninsula']
            
            # Store regional results in the dictionary
            results['IMBIE_total_mass_change_sum_east'] = IMBIE_total_mass_change_sum_east
            results['IMBIE_total_mass_change_sum_west'] = IMBIE_total_mass_change_sum_west
            results['IMBIE_total_mass_change_sum_peninsula'] = IMBIE_total_mass_change_sum_peninsula
            
            results['delta_masschange_east'] = delta_masschange_east
            results['delta_masschange_west'] = delta_masschange_west
            results['delta_masschange_peninsula'] = delta_masschange_peninsula

        # Store regional check result in the dictionary
        results['print_regionalresult_check'] = print_regionalresult_check

    return results




## Write the mass comaprision output results to csv files
def write_mass_change_comparison(icesheet, basin_result, results,mass_balance_type,start_date,end_date,csv_filename):
    print_regionalresult_check = results.get('print_regionalresult_check')

    # Data for basin mass change
    basin_mass_change_sums = basin_result['basin_mass_change_sums']
    formatted_mass_change_sums = basin_mass_change_sums.apply(lambda x: f"{x:.2f}")

    # Initialize list to store rows of data for CSV
    data_rows = []

    # Add mass change comparison header as the first row with two columns
    data_rows.append([f"Mass change comparison ({mass_balance_type})", f"{start_date} - {end_date}"])


    # Add column headers for basin mass change
    data_rows.append(['Basin', 'Model mass change (Gt)', 'IMBIE mass change (Gt)', 'Residual (Gt)'])

    # Placeholders for 'IMBIE mass change' and 'Residual' columns
    imbie_mass_change = '--'
    residual_mass_change = '--'

    # Loop through and collect each basin's subregion mass change
    for subregion, model_mass_change in formatted_mass_change_sums.items():
        data_rows.append([subregion, model_mass_change, imbie_mass_change, residual_mass_change])

    if icesheet == "AIS" and print_regionalresult_check == 'YES':
        region_mass_change_sums = basin_result.get('region_mass_change_sums')

        if region_mass_change_sums is not None:
            formatted_region_mass_change = region_mass_change_sums.apply(lambda x: f"{x:.2f}")

            # Define regions, totals, and delta changes
            regions = ["East", "West", "Peninsula", "Islands"]
            IMBIE_totals = [results.get('IMBIE_total_mass_change_sum_east'), results.get('IMBIE_total_mass_change_sum_west'),
                            results.get('IMBIE_total_mass_change_sum_peninsula')]
            delta_changes = [results.get('delta_masschange_east'), results.get('delta_masschange_west'),
                             results.get('delta_masschange_peninsula')]

            # Collect each region's mass change
            for region, total, delta in zip(regions, IMBIE_totals, delta_changes):
                mass_change = formatted_region_mass_change.get(region, "N/A")
                data_rows.append([region, mass_change, f"{total:.2f}" if total is not None else "N/A", 
                                  f"{delta:.2f}" if delta is not None else "N/A"])

    # Collect total mass balance
    IMBIE_total_mass_change_sum = results.get('IMBIE_total_mass_change_sum')
    delta_masschange = results.get('delta_masschange')
    model_total_mass_balance = basin_result['model_total_mass_balance']


    # Add the total mass balance row
    data_rows.append(['Total', f"{model_total_mass_balance:.2f}", f"{IMBIE_total_mass_change_sum:.2f}", 
                      f"{delta_masschange:.2f}"])

    # Convert the data rows into a pandas DataFrame
    df = pd.DataFrame(data_rows)

    # Write the DataFrame to a CSV file
    print(f"Writing data to CSV file: {csv_filename}")
    df.to_csv(csv_filename, index=False, header=False)



## Write the mass comaprision output results to csv files
def write_and_display_mass_change_comparison(icesheet, basin_result,results,mass_balance_type,start_date,end_date,csv_filename):
    
    print_regionalresult_check = results.get('print_regionalresult_check')
    # Data for basin mass change
    basin_mass_change_sums = basin_result['basin_mass_change_sums']
    formatted_mass_change_sums = basin_mass_change_sums.apply(lambda x: f"{x:.2f}")

    # Initialize list to store rows of data for CSV
    data_rows = []

    # Add mass change comparison header as the first row with two columns
    data_rows.append([f"Mass change comparison ({mass_balance_type})", f"{start_date} - {end_date}"])


    # Add column headers for basin mass change
    data_rows.append(['Basin', 'Model mass change (Gt)', 'IMBIE mass change (Gt)', 'Residual (Gt)'])

    # Placeholders for 'IMBIE mass change' and 'Residual' columns
    imbie_mass_change = '--'
    residual_mass_change = '--'

    print(f"\nMass change comparison ({mass_balance_type}): {start_date} - {end_date}")
    # Define column headers with fixed width for alignment
    print(f"{'Basin':<10} {'Model mass change (Gt)':<25} {'IMBIE mass change (Gt)':<25} {'Residual (Gt)':<25}")
   

    # Loop through and collect each basin's subregion mass change
    for subregion, model_mass_change in formatted_mass_change_sums.items():
        print(f"{subregion:<10} {model_mass_change:<25} {imbie_mass_change:<25} {residual_mass_change:<25}")
        data_rows.append([subregion, model_mass_change, imbie_mass_change, residual_mass_change])

    if icesheet == "AIS" and print_regionalresult_check == 'YES':
        region_mass_change_sums = basin_result.get('region_mass_change_sums')

        if region_mass_change_sums is not None:
            formatted_region_mass_change = region_mass_change_sums.apply(lambda x: f"{x:.2f}")

            # Define regions, totals, and delta changes
            regions = ["East", "West", "Peninsula", "Islands"]
            IMBIE_totals = [results.get('IMBIE_total_mass_change_sum_east'), results.get('IMBIE_total_mass_change_sum_west'),
                            results.get('IMBIE_total_mass_change_sum_peninsula')]
            delta_changes = [results.get('delta_masschange_east'), results.get('delta_masschange_west'),
                             results.get('delta_masschange_peninsula')]

            # Collect each region's mass change
            for region, total, delta in zip(regions, IMBIE_totals, delta_changes):
                mass_change = formatted_region_mass_change.get(region, "N/A")
                print(f"{region:<10} {mass_change:<25} {total:<25.2f} {delta:<25.2f}")
                data_rows.append([region, mass_change, f"{total:.2f}" if total is not None else "N/A", 
                                  f"{delta:.2f}" if delta is not None else "N/A"])

    # Collect total mass balance
    IMBIE_total_mass_change_sum = results.get('IMBIE_total_mass_change_sum')
    delta_masschange = results.get('delta_masschange')
    model_total_mass_balance = basin_result['model_total_mass_balance']


    # Print total mass balance with formatted columns
    print(f"{'Total':<10} {model_total_mass_balance.round(2):<25} {IMBIE_total_mass_change_sum:<25.2f} {delta_masschange:<25.2f}")

    
    # Add the total mass balance row
    data_rows.append(['Total', f"{model_total_mass_balance:.2f}", f"{IMBIE_total_mass_change_sum:.2f}", 
                      f"{delta_masschange:.2f}"])

    # Convert the data rows into a pandas DataFrame
    df = pd.DataFrame(data_rows)

    # Write the DataFrame to a CSV file
    print(f"Writing data to CSV file: {csv_filename}")
    df.to_csv(csv_filename, index=False, header=False)


