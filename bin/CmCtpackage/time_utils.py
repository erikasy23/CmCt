# Shared utilities
import numpy as np

def check_datarange(time_var, start_date_cftime, end_date_cftime):
    calendar_type = time_var.to_index().calendar
       
    # Get the minimum and maximum values directly from the time variable
    min_time = time_var.values.min()
    max_time = time_var.values.max()
    
    # Check if the selected start and end dates are within the range
    if min_time <= start_date_cftime <= max_time and min_time <= end_date_cftime <= max_time:
        print(f"The selected dates {{start_date_cftime}} and {{end_date_cftime}} are within the range of the model data.")
    else:
        raise ValueError(f"Error: The selected dates {{start_date_cftime}} or {{end_date_cftime}} are out of range. Model data time range is from {{min_time}} to {{max_time}}.")
