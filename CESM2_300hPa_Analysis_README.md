# CESM2 300 hPa Co-spectra Analysis

This repository contains modified MATLAB code for analyzing atmospheric momentum flux co-spectra at 300 hPa using CESM2 model output data.

## Overview

The original CESM2 analysis code has been adapted to work with the specific NetCDF files provided:
- **Historical data**: `b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc` (1961-1980)
- **RCP8.5 projection**: `b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc` (2081-2100)

## Key Modifications Made

### 1. **Data File Adaptation**
- Updated file paths to use the provided NetCDF files
- Modified variable reading to extract U and V wind components at 300 hPa level
- Adapted coordinate system handling for the specific data structure

### 2. **Pressure Level Selection**
- **300 hPa level**: Set to level 8 in the hybrid coordinate system (adjustable)
- Note: You may need to verify this level corresponds to ~300 hPa by checking actual pressure values

### 3. **Temporal Processing**
- Adapted for monthly data instead of daily data
- Implemented DJFM (December-January-February-March) season extraction
- Added interpolation/padding to create 120-day equivalent time series for spectral analysis

### 4. **Data Structure Handling**
- Added proper dimension permutation for NetCDF convention to MATLAB convention
- Implemented separate processing loops for historical and RCP8.5 data
- Updated file naming conventions for output files

### 5. **Analysis Pipeline**
- Maintained the original space-time spectral analysis methodology
- Preserved momentum flux divergence calculations
- Kept frequency-wavenumber to phase speed-wavenumber conversion
- Added improved plotting with titles and labels

## Files Included

### Main Analysis Script
- `cesm2_modified.m` - Main analysis script adapted for 300 hPa data

### Supporting Functions
- `convert_fk_to_ck.m` - Converts frequency-wavenumber to phase speed-wavenumber spectra
- `dy_dx.m` - Computes derivatives along specified dimensions
- `smooth.m` - Applies box or Gaussian smoothing windows
- `space_time_new.m` - Performs space-time spectral analysis

## Usage Instructions

### Prerequisites
- MATLAB with Signal Processing Toolbox
- NetCDF reading capabilities (`ncread` function)
- The provided NetCDF data files in the working directory

### Running the Analysis

1. **Place data files** in the MATLAB working directory:
   ```
   b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc
   b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc
   ```

2. **Run the main script**:
   ```matlab
   cesm2_modified
   ```

3. **Output files generated**:
   - `uv_DJFM_CESM2_hist_300hpa.mat` - Historical period processed data
   - `uv_DJFM_CESM2_rcp85_300hpa.mat` - RCP8.5 period processed data
   - `co-spectra_cesm2_300hpa.mat` - Final co-spectra analysis results

## Key Parameters to Verify/Adjust

### 1. **Pressure Level Selection**
```matlab
lev_300hpa = 8; % May need adjustment based on actual pressure levels
```

To verify the correct level for 300 hPa:
```matlab
% Check pressure levels (if available)
p_levels = ncread('your_file.nc', 'ps'); % or appropriate pressure variable
% Examine which level corresponds to ~300 hPa (30000 Pa)
```

### 2. **Temporal Sampling**
The code assumes monthly data and creates pseudo-daily data by replication. If your data has different temporal resolution, adjust accordingly.

### 3. **Domain and Resolution**
- **Longitude**: 360 points (1° resolution)
- **Latitude**: 180 points (1° resolution)  
- **Time**: 240 months per dataset (20 years × 12 months)

## Expected Output

The analysis produces:
1. **Co-spectra plots** showing eastward and westward propagating waves
2. **Difference plots** between RCP8.5 and historical periods
3. **Zonal wind profiles** for both time periods
4. **Phase speed-latitude diagrams** with momentum flux convergence/divergence

## Important Notes

### Data Assumptions
- **Monthly temporal resolution**: The code interpolates to daily equivalent
- **Hybrid coordinate system**: Level 8 assumed to be ~300 hPa
- **DJFM seasonal focus**: December-January-February-March analysis

### Potential Adjustments Needed
1. **Verify pressure level**: Check that level 8 corresponds to 300 hPa
2. **Temporal interpolation**: May need refinement based on actual data characteristics
3. **Domain specifics**: Adjust if your data has different spatial coverage

### Performance Considerations
- Processing time depends on data size and available RAM
- Consider reducing domain size for testing
- The script processes 20 years of data for each period

## Troubleshooting

### Common Issues
1. **File not found**: Ensure NetCDF files are in the working directory
2. **Memory errors**: Reduce domain size or process fewer years
3. **Wrong pressure level**: Verify level 8 corresponds to 300 hPa in your data

### Verification Steps
1. Check data dimensions after loading
2. Verify coordinate arrays (lat, lon, lev)
3. Examine sample wind field values for reasonableness
4. Check that seasonal extraction produces expected time series length

## Contact and Support

For questions about modifications or issues with the analysis, please check:
1. NetCDF file structure compatibility
2. MATLAB toolbox requirements
3. Memory and processing requirements for your system

---

*This analysis is adapted from atmospheric momentum flux co-spectra analysis methodology for climate model intercomparison studies.*