# CESM2 Co-Spectra Analysis at 300 hPa

## Overview
This code performs space-time spectral analysis of atmospheric momentum flux divergence using CESM2 model output data. The analysis focuses on eastward and westward propagating waves at the 300 hPa pressure level.

## Data Files Required
- `b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc` (Historical period: 1961-1980)
- `b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc` (Future period: 2081-2100 under RCP8.5)

## Key Modifications Made

### 1. Data Reading and Level Selection
- Modified to read U and V wind components at 300 hPa (level index 16)
- Updated file paths to match the provided NetCDF files
- Reads data using `ncread` with proper indexing: `[lon, lev, lat, time]`

### 2. Temporal Coverage
- **Original**: Used winter seasons only (DJFM, 120 days)
- **Modified**: Uses all monthly data (12 months per year × 20 years = 240 months)
- Applied Hanning window with 12-month length instead of 120-day length

### 3. Data Processing
- Processes both historical (1961-1980) and RCP8.5 (2081-2100) datasets
- Maintains the same space-time spectral analysis methodology
- Calculates momentum flux divergence and applies smoothing

### 4. Output
- Generates phase speed vs. latitude plots showing differences between future and historical periods
- Saves intermediate results in `.mat` files for later analysis
- Creates contour plots with proper color mapping

## Files Structure

### Main Script
- `cesm2_cospectra_analysis.m` - Main analysis script

### Supporting Functions
- `convert_fk_to_ck.m` - Converts frequency-wavenumber to phase speed-wavenumber space
- `dy_dx.m` - Calculates derivatives along specified dimensions
- `smooth.m` - Provides smoothing functions (box and Gaussian)
- `space_time_new.m` - Performs space-time spectral analysis

## Running the Code

1. Ensure all NetCDF files are in the `/workspace/` directory
2. Ensure all MATLAB functions are in the same directory
3. Run the main script: `cesm2_cospectra_analysis.m`
4. The code will:
   - Read and process the data
   - Perform space-time analysis
   - Generate and save results
   - Create visualization plots

## Key Parameters
- **Pressure Level**: 300 hPa (index 16 in the level array)
- **Time Window**: 12 months (full year)
- **Spatial Coverage**: Global (360 longitude × 180 latitude)
- **Temporal Coverage**: 20 years for each period

## Output Files
- `uv_monthly_CESM2_hist_300hpa.mat` - Historical period results
- `uv_monthly_CESM2_rcp85_300hpa.mat` - RCP8.5 period results
- `co-spectra_cesm2_300hpa.mat` - Final co-spectra results

## Physical Interpretation
The analysis computes the co-spectra of momentum flux divergence, which helps understand:
- Wave propagation characteristics in the atmosphere
- Changes in atmospheric dynamics between historical and future climate
- Eastward vs. westward propagating wave patterns
- Latitudinal distribution of momentum flux changes

## Requirements
- MATLAB with Signal Processing Toolbox (for `hann` function)
- Access to `brewermap` function for color mapping (or substitute with built-in colormaps)
- Sufficient memory for processing large 4D arrays
