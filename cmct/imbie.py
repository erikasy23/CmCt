from .time_utils import check_datarange,days_in_year
import os
import numpy as np
import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import datetime
from datetime import timedelta 




### Load the model data and calculate model mass balance for each basin and total mass balance for whole region
## Interpolate the data to each IMBIE time and calculate the time varying mass change
def process_model_data(mod_ds,time_var, IMBIE_total_mass_change_sum, \
                       start_date_cftime, end_date_cftime, start_date_fract, end_date_fract, \
                       rho_ice, projection, shape_filename, icesheet):
    
    # Model data
    lithk = mod_ds['lithk']
    
    lithk['time'] = [date.year + (date.dayofyr-1) / days_in_year(date) for date in time_var.data]

    # Load basin shapefile 
    basins_gdf = gpd.read_file(shape_filename)
    
    # Check the selcted dates are within the range of model data
    check_datarange(time_var,start_date_cftime, end_date_cftime)
        
    # Set start_date as the first date in 'Year' and filtered_time_var as all subsequent dates
    start_date_imbie = start_date_fract
    filtered_time_var = IMBIE_total_mass_change_sum['Year'].iloc[0:]
    
    #calculate area = x_resolution*y_resolution
    x_coords = mod_ds['x'].values
    y_coords = mod_ds['y'].values
    x_resolution = abs(x_coords[1] - x_coords[0])
    y_resolution = abs(y_coords[1] - y_coords[0])
    
    # Create a list of Point geometries from coordinate grids
    points = [Point(x, y) for x in x_coords for y in y_coords]
    
    # Initialize a dictionary to store residuals
    model_mass_change = {}
    
    # Interpolate limnsw at the start date (initial reference)
    lithk_start = lithk.interp(time=start_date_fract).data.transpose().flatten()
    
    # Loop through each filtered time step to calculate the residual
    for i, time_step in enumerate(filtered_time_var):
     
        # Interpolate limnsw at the current time_step
        lithk_current = lithk.interp(time=time_step).data.transpose().flatten()
        
        # Calculate the residual (difference from start)
        lithk_delta = lithk_current - lithk_start
    
        lithk_delta[np.isnan(lithk_delta)] = 0
    
        #calculate area = x_resolution*y_resolution
        lithk_delta = (lithk_delta * x_resolution*y_resolution)*rho_ice * 1e-12

        ## TOTAL AREA
        # Sum all of the area mass change
        model_total_mass_balance_unmasked= np.nansum(lithk_delta)
               
        ## BASIN AREA
        # Sum all of the basin mass change
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
        model_total_mass_balance_masked = basin_mass_change_sums.sum()
          

        # Store the residual in the dictionary for the current time step
        model_mass_change[str(time_step)] = {
            'model_total_mass_balance_unmasked': model_total_mass_balance_unmasked,
            'model_total_mass_balance_masked': model_total_mass_balance_masked,
            'basin_mass_change_sums': basin_mass_change_sums,
            'region_mass_change_sums': region_mass_change_sums
        } 
           
        
    # Return all results as a dictionary
    return model_mass_change




### Extract time varying IMBIE mass balance data and calculate the time varying mass difference 
def process_imbie_data(obs_filename,start_date_fract,end_date_fract,mass_balance_column):

    # Load the CSV file
    mass_balance_data = pd.read_csv(obs_filename)
    
    # Column names
    date_column = 'Year'
    
    # Ensure the 'Year' column is treated as float to capture the fractional year part
    mass_balance_data['Year'] = mass_balance_data['Year'].astype(float)
    

    # Sort the data by 'Date' column to ensure it’s in increasing order of both year and fraction
    mass_balance_data = mass_balance_data.sort_values(by='Year')
        

    
    # Check if the column exists in the DataFrame
    if mass_balance_column not in mass_balance_data.columns:
        raise ValueError(f"Error: The column '{mass_balance_column}' does not exist in the CSV file.")
    
    
    # Get the initial mass balance value for the start date
    data_start_date = mass_balance_data[mass_balance_data['Year'] == start_date_fract]
    if data_start_date.empty:
        raise ValueError(f"Error: No data available for the start date {start_date_fract}.")
    mass_balance_start_value = data_start_date[mass_balance_column].iloc[0]  # value of start date
    
    # Filter data between start_date_converted and end_date_converted (inclusive)
    filtered_data = mass_balance_data[
        (mass_balance_data['Year'] > start_date_fract) & (mass_balance_data['Year'] <= end_date_fract)
    ].copy()
    
    
    # Initialize the previous day's mass balance value to the starting mass balance
    previous_mass_balance = mass_balance_start_value
    
    # Calculate daily mass change from the previous day for each time step
    mass_changes = []  # To store the daily mass changes
    
    for index, row in filtered_data.iterrows():
        current_mass_balance = row[mass_balance_column]
        # Calculate the change from the previous day's balance
        mass_change = current_mass_balance-previous_mass_balance
        mass_changes.append(mass_change)
   
    
    # Assign the calculated mass changes to a new column in the DataFrame
    filtered_data['IMBIE_Mass_Change'] = mass_changes


    return filtered_data
   
    
 

### Calculate mass balance difference of IMBIE and model data
def calculate_model_imbie_residuals(start_date_fract, end_date_fract, \
                  icesheet, basin_result, IMBIE_total_mass_change_sum, mass_balance_column, \
                  obs_east_filename=None, obs_west_filename=None, obs_peninsula_filename=None):
    # Initialize a dictionary to store results
    results = {}
    print_regionalresult_check = 'NO'
    
    # Loop through IMBIE_total_mass_change_sum to populate results
    for i, row in IMBIE_total_mass_change_sum.iterrows():
        date = row['Year']  # Ensure this matches the date format in basin_result keys
        imbie_mass_change = row['IMBIE_Mass_Change']
        
        # Check if the date exists in basin_result
        if str(date) in basin_result:
            model_mass_change_masked = basin_result[str(date)]['model_total_mass_balance_masked']
            model_mass_change_unmasked = basin_result[str(date)]['model_total_mass_balance_unmasked']
            # Calculate the delta
            delta_masschange_masked  = imbie_mass_change - model_mass_change_masked 
            delta_masschange_unmasked  = imbie_mass_change - model_mass_change_unmasked 
            
            # Store results in the dictionary
            results[str(date)] = {
                'IMBIE_total_mass_change_sum': imbie_mass_change,
                'Delta_MassChange_masked': delta_masschange_masked,
                'Delta_MassChange_unmasked': delta_masschange_unmasked
            }
        else:
            print(f"Date {date} not found in basin_result. Skipping.")  
    

    if icesheet == "AIS":  
        # Check if the required files exist
        if (obs_east_filename and os.path.exists(obs_east_filename)) and \
           (obs_west_filename and os.path.exists(obs_west_filename)) and \
           (obs_peninsula_filename and os.path.exists(obs_peninsula_filename)):

            print_regionalresult_check = 'YES'
            
            # Calculate total mass for each region
            IMBIE_total_mass_change_sum_east = process_imbie_data(obs_east_filename, start_date_fract,end_date_fract, mass_balance_column)
            IMBIE_total_mass_change_sum_west = process_imbie_data(obs_west_filename, start_date_fract,end_date_fract, mass_balance_column)
            IMBIE_total_mass_change_sum_peninsula = process_imbie_data(obs_peninsula_filename, start_date_fract,end_date_fract, mass_balance_column)
    
            # Loop through IMBIE_total_mass_change_sum to populate results for each region
            regional_results = {}

            # Assuming the datasets have consistent row structures and indices
            for i, (east_row, west_row, peninsula_row) in enumerate(zip(
                IMBIE_total_mass_change_sum_east.iterrows(),
                IMBIE_total_mass_change_sum_west.iterrows(),
                IMBIE_total_mass_change_sum_peninsula.iterrows()
            )):
                date = east_row[1]['Year']  # Assuming the 'Year' column is consistent across datasets
            
                # Extract values for each region
                imbie_mass_change_east = east_row[1]['IMBIE_Mass_Change']
                imbie_mass_change_west = west_row[1]['IMBIE_Mass_Change']
                imbie_mass_change_peninsula = peninsula_row[1]['IMBIE_Mass_Change']
                
                
                # Check if the date exists in basin_result
                if str(date) in basin_result:
                    region_mass_change_sums = basin_result[str(date)]['region_mass_change_sums']
                    
                    # East region
                    if 'East' in region_mass_change_sums:
                        delta_masschange_east = imbie_mass_change_east - region_mass_change_sums['East']
                        regional_results.setdefault(str(date), {}).update({
                            'IMBIE_Mass_Change_East': imbie_mass_change_east,
                            'Delta_MassChange_East': delta_masschange_east
                        })
                    
                    # West region
                    if 'West' in region_mass_change_sums:
                        delta_masschange_west = imbie_mass_change_west - region_mass_change_sums['West']
                        regional_results.setdefault(str(date), {}).update({
                            'IMBIE_Mass_Change_West': imbie_mass_change_west,
                            'Delta_MassChange_West': delta_masschange_west
                        })
                    
                    # Peninsula region
                    if 'Peninsula' in region_mass_change_sums:
                        delta_masschange_peninsula = imbie_mass_change_peninsula - region_mass_change_sums['Peninsula']
                        regional_results.setdefault(str(date), {}).update({
                            'IMBIE_Mass_Change_Peninsula': imbie_mass_change_peninsula,
                            'Delta_MassChange_Peninsula': delta_masschange_peninsula
                        })
                else:
                    print(f"Date {date} not found in basin_result. Skipping.")
            
            # Store regional results in the main results dictionary
            results['Regional_Mass_Change_Summary'] = regional_results

    # Store regional check result in the dictionary
    results['print_regionalresult_check'] = print_regionalresult_check
    return results




def write_and_display_mass_change_comparison_all_dates(icesheet, basin_result, results, mass_balance_type,start_date_fract,end_date_fract, csv_filename):
    # Initialize list to store rows of data for CSV
    data_rows = []

    print_regionalresult_check = results.get('print_regionalresult_check')
    
    # Add mass change comparison header
    data_rows.append([f"Mass change comparison ({mass_balance_type})", f"{start_date_fract} - {end_date_fract}"])
    data_rows.append(['Date', 'Basin/Region', 'Model mass change (Gt)', 'IMBIE mass change (Gt)', 'Residual (Gt)'])

    # Determine basins and regions from the first available date after the start_date
    basins = []
    regions = []
    for date in sorted(basin_result.keys()):
        if float(date) > start_date_fract:
            basins = list(basin_result[date].get('basin_mass_change_sums', {}).keys())
            if icesheet == "AIS" and print_regionalresult_check == 'YES':
                regions = list(basin_result[date].get('region_mass_change_sums', {}).keys())
            break

    # Add rows for each basin with zero values for the start_date
    for basin in basins:
        data_rows.append([start_date_fract, basin, "0.00", "--", "--"])

    # Add rows for each region with zero values for the start_date (if applicable)
    if icesheet == "AIS" and print_regionalresult_check == 'YES':
        for region in regions:
            data_rows.append([start_date_fract, region, "0.00", "0.00", "0.00"])

    # Add totals (masked and unmasked) with zero values for the start_date
    data_rows.append([start_date_fract, 'Masked_Total', "0.00", "0.00", "0.00"])
    data_rows.append([start_date_fract, 'Unmasked_Total', "0.00", "0.00", "0.00"])
 

    print(f"\n Time-varying Mass change comparison ({mass_balance_type}): {start_date_fract} - {end_date_fract}")
    print(f"{'Date':<15} {'Basin/Region':<20} {'Model mass change (Gt)':<25} {'IMBIE mass change (Gt)':<25} {'Residual (Gt)':<20}")
    model_total_mass_balance_masked=0.00
    imbie_total_mass_change_sum=0.00
    delta_masschange_masked=0.00
    
    print(f"{start_date_fract:<15} {'Masked_Total':<20} {model_total_mass_balance_masked :<25} {imbie_total_mass_change_sum:<25} {delta_masschange_masked :<20}")
    

    # Process the rest of the dates
    for date, result in results.items():
        if date in basin_result:
            # Basin mass change sums
            basin_mass_change_sums = basin_result[date].get('basin_mass_change_sums', {})

            for basin, model_mass_change in basin_mass_change_sums.items():
                imbie_mass_change = '--'
                residual_mass_change = '--'
                data_rows.append([date, basin, f"{model_mass_change:.2f}", imbie_mass_change, residual_mass_change])

            if icesheet == "AIS" and print_regionalresult_check == 'YES':
                # Regional mass change sums
                region_mass_change_sums = basin_result[date].get('region_mass_change_sums', {})
                for region, model_mass_change in region_mass_change_sums.items():
                    imbie_mass_change = results.get('Regional_Mass_Change_Summary', {}).get(date, {}).get(f'IMBIE_Mass_Change_{region}', '--')
                    residual_mass_change = results.get('Regional_Mass_Change_Summary', {}).get(date, {}).get(f'Delta_MassChange_{region}', '--')
        
                    imbie_mass_change = f"{imbie_mass_change:.2f}" if isinstance(imbie_mass_change, (float, int)) else "--"
                    residual_mass_change = f"{residual_mass_change:.2f}" if isinstance(residual_mass_change, (float, int)) else "--"
                    data_rows.append([date, region, f"{model_mass_change:.2f}", imbie_mass_change, residual_mass_change])
            
            # Total mass balance masked
            model_total_mass_balance_masked = basin_result[date].get('model_total_mass_balance_masked', '--')   
            imbie_total_mass_change_sum = result.get('IMBIE_total_mass_change_sum', '--')
            delta_masschange_masked  = result.get('Delta_MassChange_masked', '--')
    
            model_total_mass_balance_masked  = f"{model_total_mass_balance_masked :.2f}" if isinstance(model_total_mass_balance_masked , (float, int)) else "--"
            imbie_total_mass_change_sum = f"{imbie_total_mass_change_sum:.2f}" if isinstance(imbie_total_mass_change_sum, (float, int)) else "--"
            delta_masschange_masked  = f"{delta_masschange_masked :.2f}" if isinstance(delta_masschange_masked , (float, int)) else "--"
            data_rows.append([date, 'Masked_Total', model_total_mass_balance_masked , imbie_total_mass_change_sum, delta_masschange_masked ]) 
            print(f"{date:<15} {'Masked_Total':<20} {model_total_mass_balance_masked :<25} {imbie_total_mass_change_sum:<25} {delta_masschange_masked :<20}")

            # Total mass balance unmasked
            model_total_mass_balance_unmasked = basin_result[date].get('model_total_mass_balance_unmasked', '--')
            imbie_total_mass_change_sum = result.get('IMBIE_total_mass_change_sum', '--')
            delta_masschange_unmasked = result.get('Delta_MassChange_unmasked', '--')
        
        
            model_total_mass_balance_unmasked = f"{model_total_mass_balance_unmasked:.2f}" if isinstance(model_total_mass_balance_unmasked, (float, int)) else "--"
            imbie_total_mass_change_sum = f"{imbie_total_mass_change_sum:.2f}" if isinstance(imbie_total_mass_change_sum, (float, int)) else "--"
            delta_masschange_unmasked = f"{delta_masschange_unmasked:.2f}" if isinstance(delta_masschange_unmasked, (float, int)) else "--"
            data_rows.append([date, 'Unmasked_Total', model_total_mass_balance_unmasked, imbie_total_mass_change_sum, delta_masschange_unmasked])

    # Convert the data rows into a pandas DataFrame
    df = pd.DataFrame(data_rows)
    
    # Write the DataFrame to a CSV file
    print(f"\nWriting data to CSV file: {csv_filename}")
    df.to_csv(csv_filename, index=False, header=False)




def write_mass_change_comparison_all_dates(icesheet, basin_result, results, mass_balance_type,start_date_fract,end_date_fract, csv_filename):
    # Initialize list to store rows of data for CSV
    data_rows = []

    print_regionalresult_check = results.get('print_regionalresult_check')
    
    # Add mass change comparison header
    data_rows.append([f"Mass change comparison ({mass_balance_type})", f"{start_date_fract} - {end_date_fract}"])
    data_rows.append(['Date', 'Basin/Region', 'Model mass change (Gt)', 'IMBIE mass change (Gt)', 'Residual (Gt)'])

    # Determine basins and regions from the first available date after the start_date
    basins = []
    regions = []
    for date in sorted(basin_result.keys()):
        if float(date) > start_date_fract:
            basins = list(basin_result[date].get('basin_mass_change_sums', {}).keys())
            if icesheet == "AIS" and print_regionalresult_check == 'YES':
                regions = list(basin_result[date].get('region_mass_change_sums', {}).keys())
            break

    # Add rows for each basin with zero values for the start_date
    for basin in basins:
        data_rows.append([start_date_fract, basin, "0.00", "--", "--"])

    # Add rows for each region with zero values for the start_date (if applicable)
    if icesheet == "AIS" and print_regionalresult_check == 'YES':
        for region in regions:
            data_rows.append([start_date_fract, region, "0.00", "0.00", "0.00"])

    # Add totals (masked and unmasked) with zero values for the start_date
    data_rows.append([start_date_fract, 'Masked_Total', "0.00", "0.00", "0.00"])
    data_rows.append([start_date_fract, 'Unmasked_Total', "0.00", "0.00", "0.00"])
 

    model_total_mass_balance_masked=0.00
    imbie_total_mass_change_sum=0.00
    delta_masschange_masked=0.00
    

    # Process the rest of the dates
    for date, result in results.items():
        if date in basin_result:
            # Basin mass change sums
            basin_mass_change_sums = basin_result[date].get('basin_mass_change_sums', {})

            for basin, model_mass_change in basin_mass_change_sums.items():
                imbie_mass_change = '--'
                residual_mass_change = '--'
                data_rows.append([date, basin, f"{model_mass_change:.2f}", imbie_mass_change, residual_mass_change])

            if icesheet == "AIS" and print_regionalresult_check == 'YES':
                # Regional mass change sums
                region_mass_change_sums = basin_result[date].get('region_mass_change_sums', {})
                for region, model_mass_change in region_mass_change_sums.items():
                    imbie_mass_change = results.get('Regional_Mass_Change_Summary', {}).get(date, {}).get(f'IMBIE_Mass_Change_{region}', '--')
                    residual_mass_change = results.get('Regional_Mass_Change_Summary', {}).get(date, {}).get(f'Delta_MassChange_{region}', '--')
        
                    imbie_mass_change = f"{imbie_mass_change:.2f}" if isinstance(imbie_mass_change, (float, int)) else "--"
                    residual_mass_change = f"{residual_mass_change:.2f}" if isinstance(residual_mass_change, (float, int)) else "--"
                    data_rows.append([date, region, f"{model_mass_change:.2f}", imbie_mass_change, residual_mass_change])
            
            # Total mass balance masked
            model_total_mass_balance_masked = basin_result[date].get('model_total_mass_balance_masked', '--')   
            imbie_total_mass_change_sum = result.get('IMBIE_total_mass_change_sum', '--')
            delta_masschange_masked  = result.get('Delta_MassChange_masked', '--')
    
            model_total_mass_balance_masked  = f"{model_total_mass_balance_masked :.2f}" if isinstance(model_total_mass_balance_masked , (float, int)) else "--"
            imbie_total_mass_change_sum = f"{imbie_total_mass_change_sum:.2f}" if isinstance(imbie_total_mass_change_sum, (float, int)) else "--"
            delta_masschange_masked  = f"{delta_masschange_masked :.2f}" if isinstance(delta_masschange_masked , (float, int)) else "--"
            data_rows.append([date, 'Masked_Total', model_total_mass_balance_masked , imbie_total_mass_change_sum, delta_masschange_masked ]) 

            # Total mass balance unmasked
            model_total_mass_balance_unmasked = basin_result[date].get('model_total_mass_balance_unmasked', '--')
            imbie_total_mass_change_sum = result.get('IMBIE_total_mass_change_sum', '--')
            delta_masschange_unmasked = result.get('Delta_MassChange_unmasked', '--')
        
        
            model_total_mass_balance_unmasked = f"{model_total_mass_balance_unmasked:.2f}" if isinstance(model_total_mass_balance_unmasked, (float, int)) else "--"
            imbie_total_mass_change_sum = f"{imbie_total_mass_change_sum:.2f}" if isinstance(imbie_total_mass_change_sum, (float, int)) else "--"
            delta_masschange_unmasked = f"{delta_masschange_unmasked:.2f}" if isinstance(delta_masschange_unmasked, (float, int)) else "--"
            data_rows.append([date, 'Unmasked_Total', model_total_mass_balance_unmasked, imbie_total_mass_change_sum, delta_masschange_unmasked])

    # Convert the data rows into a pandas DataFrame
    df = pd.DataFrame(data_rows)
    
    # Write the DataFrame to a CSV file
    print(f"\nWriting data to CSV file: {csv_filename}")
    df.to_csv(csv_filename, index=False, header=False)




