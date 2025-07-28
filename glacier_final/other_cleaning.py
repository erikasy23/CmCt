import os
import pandas as pd
import xarray as xr
import numpy as np
import io

# ======================================================
# ASSUMPTIONS & DATA REQUIREMENTS FOR OTHER_CLEANING
# ======================================================
#
# 1. GENERAL
# - The user must define:
#     - source_type: Either "Upload" (CSV) or "Link" (.nc files).
#     - selected_model: The model name to process (e.g., "OGGM", "custom_model").
#     - selected_scenario: The scenario name to process (e.g., "historical").
#     - selected_regions: A list of region codes (e.g., ["RGI01", "RGI06"]).
# - The script processes only data matching (selected_model, selected_scenario).
# - Output: A DataFrame (summary_df) with columns:
#       ['calendar_year', 'dmdt', 'dhdt',
#        'dmdt_5yr_avg', 'dhdt_5yr_avg',
#        'model', 'scenario', 'region'].
#
# ------------------------------------------------------
# 2. CSV UPLOADS (source_type == "Upload")
# - File type: A single .csv file.
# - Required columns:
#       'model'        -> Model name (must match selected_model).
#       'scenario'     -> Scenario name (must match selected_scenario).
#       'region'       -> Region code for each glacier entry.
#       'rgi_id'       -> Unique glacier identifiers.
#       'volume'       -> Glacier volume values (for dmdt).
#       'area'         -> Glacier area values (for dhdt).
#       'calendar_year' OR 'hydro_year' -> One must be present.
# - Auto-conversion:
#       If only 'hydro_year' exists, 'calendar_year' is generated (hydro_year - 1).
#       If only 'calendar_year' exists, 'hydro_year' is generated (calendar_year + 1).
# - Filtering:
#       Rows are filtered by (model == selected_model) & (scenario == selected_scenario).
#       Data is then looped over selected_regions for further processing.
#
# ------------------------------------------------------
# 3. NETCDF FILES (source_type == "Link")
# - File type: A directory containing multiple .nc files.
# - File naming format:
#       Filenames must contain the model and scenario as parts [5] and [6] when split by "_".
#       Example: some_output_run_v1_modelA_historical_RGI01.nc
#                parts[5] = "modelA", parts[6] = "historical".
# - Variables inside .nc files:
#       'volume'  -> Glacier volume (time × rgi_id).
#       'area'    -> Glacier area (time × rgi_id).
#       'time'    -> Represents calendar_year.
# - Filtering:
#       Only files matching (selected_model, selected_scenario) are processed.
#       Regions are identified by checking if region string (e.g., "RGI01") is in the filename.
#
# ------------------------------------------------------
# 4. CALCULATIONS
# - Mass change (dmdt):
#       dmdt = Δvolume * 900 (kg/m³) / 1e12  -> Gt/year.
# - Elevation change (dhdt):
#       dhdt = Δvolume / area                -> m/year.
# - 5-Year Averages:
#       Data is grouped into 5-year blocks.
#       The average dmdt and dhdt are assigned to the last year of each block.
#       All other years within the block are set to NaN for these averages.
#
# ------------------------------------------------------
# 5. OUTPUT
# - The final summary_df contains:
#       'calendar_year', 'dmdt', 'dhdt',
#       'dmdt_5yr_avg', 'dhdt_5yr_avg',
#       'model', 'scenario', 'region'.
# - This summary_df is used downstream for plotting and comparisons.
#
# ======================================================


# ======================================================
# Helper Functions
# ======================================================

def validate_required_columns(df, cols):
    """Ensure DataFrame has required columns."""
    missing = cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

def ensure_calendar_and_hydro_year(df):
    """Ensure both calendar_year and hydro_year are present."""
    if 'calendar_year' not in df.columns and 'hydro_year' not in df.columns:
        raise ValueError("Data must include either 'calendar_year' or 'hydro_year'.")
    if 'calendar_year' not in df.columns:
        df['calendar_year'] = df['hydro_year'] - 1
    if 'hydro_year' not in df.columns:
        df['hydro_year'] = df['calendar_year'] + 1
    return df

def load_csv_from_upload(uploaded_file):
    """Load CSV data, validate, and filter by selected model/scenario."""
    content = io.BytesIO(uploaded_file["content"])
    df = pd.read_csv(content)

    required_cols = {"model", "scenario", "region", "rgi_id", "volume", "area"}
    validate_required_columns(df, required_cols)
    df = ensure_calendar_and_hydro_year(df)

    # Filter by selected model and scenario
    df = df[(df["model"] == selected_model) & (df["scenario"] == selected_scenario)]
    if df.empty:
        raise ValueError(f"No data found for model={selected_model}, scenario={selected_scenario} in CSV.")
    
    return df

def collect_model_scenario_files(link_path):
    """
    Collect files by (model, scenario).
    Returns a dict: { (model, scenario): [file1, file2, ...] }
    """
    all_files = [f for f in os.listdir(link_path) if f.endswith(".nc")]
    if not all_files:
        raise ValueError(f"No .nc files found in {link_path}")

    model_scenario_map = {}
    for fname in all_files:
        parts = fname.split("_")
        if len(parts) < 7:
            print(f"Skipping {fname}: does not follow OGGM naming convention.")
            continue
        model = parts[5]
        scenario = parts[6]
        model_scenario_map.setdefault((model, scenario), []).append(fname)
    return model_scenario_map

def compute_mass_elevation_changes(ds):
    """
    Compute dmdt and dhdt from xarray Dataset.
    """
    volume = ds["volume"]
    area = ds["area"]

    # Mass change (Gt/year)
    mass = volume * 900  # kg
    mass_change = mass.diff(dim="time") / 1e12  # Gt/year

    # Elevation change (m/year)
    volume_change = volume.diff(dim="time")
    elevation_change = volume_change / area

    # Handle invalid cases
    invalid = ((area == 0) & (volume_change != 0)) | ((area != 0) & (volume_change == 0))
    elevation_change = elevation_change.where(~invalid, np.nan)
    both_zero = (area == 0) & (volume_change == 0)
    elevation_change = elevation_change.where(~both_zero, 0)

    return elevation_change, mass_change

def compute_5yr_averages(df):
    """Add 5-year averaged columns for dmdt and dhdt."""
    block_size = 5
    calendar_blocks = (df['calendar_year'] - df['calendar_year'].min()) // block_size
    df['dmdt_5yr_avg'] = df.groupby(calendar_blocks)['dmdt'].transform('mean')
    df['dhdt_5yr_avg'] = df.groupby(calendar_blocks)['dhdt'].transform('mean')

    # Set all but last year of block to NaN
    for block in df.groupby(calendar_blocks).groups.values():
        df.loc[block[:-1], ['dmdt_5yr_avg', 'dhdt_5yr_avg']] = np.nan

    return df

def process_nc_combination(link_path, files, model, scenario, regions):
    """
    Load and process all .nc files for a single (model, scenario) combination.
    Loops through regions, producing summary_df for each.
    """
    summary_frames = []
    for region in regions:
        ds_list = [xr.open_dataset(os.path.join(link_path, f)) for f in files if region in f]
        if not ds_list:
            print(f"No data found for region={region} in model={model}, scenario={scenario}")
            continue

        combined_ds = xr.concat(ds_list, dim="rgi_id")

        # Compute changes
        elevation_change, mass_change = compute_mass_elevation_changes(combined_ds)

        # Aggregate by year
        dhdt_simple_avg = elevation_change.mean(dim="rgi_id")
        dmdt_sum = mass_change.sum(dim="rgi_id")

        # Build summary DataFrame
        summary_df = pd.DataFrame({
            "calendar_year": combined_ds["time"].values[1:],  # after diff
            "dmdt": dmdt_sum.values,
            "dhdt": dhdt_simple_avg.values,
            "model": model,
            "scenario": scenario,
            "region": region
        })
        summary_df = compute_5yr_averages(summary_df)
        summary_frames.append(summary_df)

    return pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame()

# ======================================================
# Main Cleaning Logic
# ======================================================

custom_model_data = None
model_to_scenarios = {}
summary_df = pd.DataFrame()

if source_type == "Upload":
    custom_model_data = load_csv_from_upload(uploaded_file)

    # Loop through selected regions
    summary_frames = []
    for region in selected_regions:
        df_region = custom_model_data[custom_model_data["region"] == region]
        if df_region.empty:
            print(f"No data for region {region} in CSV.")
            continue

        # Calculate mass and elevation change
        df_region = df_region.sort_values(["region", "calendar_year", "rgi_id"])
        df_region["dmdt"] = df_region.groupby("rgi_id")["volume"].diff() * 900 / 1e12
        df_region["dhdt"] = df_region.groupby("rgi_id")["volume"].diff() / df_region["area"]

        df_region = compute_5yr_averages(df_region)
        summary_frames.append(df_region)

    summary_df = pd.concat(summary_frames, ignore_index=True)

elif source_type == "Link":
    model_scenario_map = collect_model_scenario_files(link_path)
    # Filter to selected model and scenario
    if (selected_model, selected_scenario) not in model_scenario_map:
        raise ValueError(f"No files found for model={selected_model}, scenario={selected_scenario}")
    files = model_scenario_map[(selected_model, selected_scenario)]
    summary_df = process_nc_combination(link_path, files, selected_model, selected_scenario, selected_regions)

else:
    raise ValueError("Invalid source_type for 'Other (Custom)' model. Must be 'Upload' or 'Link'.")

print("Custom model data loaded and cleaned successfully.")
print(f"Processed model={selected_model}, scenario={selected_scenario}")
# Save the cleaned data to a CSV file
output_file_path = 'cleaned_custom_model_data.csv'
os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
summary_df.to_csv(output_file_path, index=False)