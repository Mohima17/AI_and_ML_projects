# CESM2 Co-Spectra Analysis at 300 hPa

This code performs space-time spectral analysis of atmospheric wind data from CESM2 NetCDF files at the 300 hPa pressure level.

## Files Included

### Main Analysis Script
- `cesm2_300hpa_analysis.m` - Modified main script for analyzing NetCDF data at 300 hPa

### Supporting Functions
- `convert_fk_to_ck.m` - Converts frequency-wavenumber to phase speed-wavenumber spectra
- `dy_dx.m` - Computes derivatives along specified dimensions (2D and 3D)
- `smooth.m` - Smoothing function with box and Gaussian window options
- `space_time_new.m` - Space-time spectral analysis function

## Key Modifications Made

1. **NetCDF Data Reading**: 
   - Modified to read U and V wind components from your specific NetCDF files
   - Uses level index 16 for 300 hPa as specified
   - Reads both historical (1961-1980) and RCP8.5 (2081-2100) data

2. **Data Structure Adaptation**:
   - Adjusted for your data's temporal resolution (monthly data)
   - Modified DJFM (December-January-February-March) extraction logic
   - Handles 20-year periods for both historical and future scenarios

3. **Variable Naming**:
   - Updated variable names to distinguish between historical and RCP data
   - Added `_300hpa` suffix to output file names for clarity

4. **Flexible Array Sizing**:
   - Made wavenumber and frequency arrays adaptive to actual data dimensions
   - Added safety checks for array bounds

## Input Data Requirements

Your NetCDF files should contain:
- **Variables**: U (zonal wind), V (meridional wind), lat, lon, lev, time
- **Dimensions**: time, lev, lat, lon
- **Level 16**: Should correspond to 300 hPa pressure level
- **Time coverage**: Monthly data for the specified periods

### File Names Expected:
- `b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc` (Historical)
- `b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc` (RCP8.5)

## Usage

1. Ensure your NetCDF files are in the working directory
2. Run the main script:
   ```matlab
   cesm2_300hpa_analysis
   ```

## Output

The script produces:
- Intermediate data files: `uv_DJFM_CESM2_hist_300hpa.mat` and `uv_DJFM_CESM2_rcp85_300hpa.mat`
- Final results: `co-spectra_cesm2_300hpa.mat`
- Visualization: Co-spectra difference plot showing changes between historical and future periods

## Dependencies

- MATLAB with NetCDF reading capabilities
- Signal Processing Toolbox (for `hann` function)
- Statistics and Machine Learning Toolbox (if using `brewermap` colormap function)

## Notes

- The code assumes monthly data and extracts DJFM months for analysis
- Level index 16 is hardcoded as the 300 hPa level
- The supporting MATLAB scripts (`convert_fk_to_ck.m`, `dy_dx.m`, `smooth.m`, `space_time_new.m`) do not require modifications and work as originally designed
- Adjust the `nyears_hist` and `nyears_rcp` variables if your data covers different time periods
