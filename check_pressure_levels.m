% Helper script to check pressure levels and find 300 hPa level
% Run this before the main analysis to verify the correct level index

clear; close all;

% File names
hist_file = 'b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc';

% Read coordinate variables
fprintf('Reading coordinate information...\n');
lat = double(ncread(hist_file, 'lat'));
lon = double(ncread(hist_file, 'lon'));
lev = double(ncread(hist_file, 'lev'));
time = double(ncread(hist_file, 'time'));

fprintf('Data dimensions:\n');
fprintf('  Longitude: %d points (%.1f to %.1f)\n', length(lon), min(lon), max(lon));
fprintf('  Latitude: %d points (%.1f to %.1f)\n', length(lat), min(lat), max(lat));
fprintf('  Levels: %d levels\n', length(lev));
fprintf('  Time: %d time steps\n', length(time));

% Display level information
fprintf('\nLevel information:\n');
fprintf('Level Index | Level Value\n');
fprintf('------------|------------\n');
for i = 1:length(lev)
    fprintf('    %2d      |   %.6f\n', i, lev(i));
end

% Try to read pressure field to determine actual pressure levels
try
    fprintf('\nReading surface pressure field...\n');
    ps = ncread(hist_file, 'PS', [1, 1, 1], [Inf, Inf, 1]); % Read first time step
    ps_sample = squeeze(ps(180, 90)); % Sample point (near equator, middle longitude)
    
    fprintf('Sample surface pressure: %.2f Pa (%.2f hPa)\n', ps_sample, ps_sample/100);
    
    % For hybrid coordinates, pressure at level k is approximately:
    % p(k) = a(k) + b(k) * ps
    % Since we don't have a and b coefficients, we'll estimate
    
    fprintf('\nEstimated pressure levels (assuming typical atmospheric profile):\n');
    fprintf('Level | Estimated Pressure (hPa)\n');
    fprintf('------|----------------------\n');
    
    % Typical pressure levels for CAM model (approximate)
    typical_pressures = [3, 7, 14, 25, 39, 61, 89, 127, 177, 239, 319, 424, 551, 704, 867, ...
                        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000];
    
    for i = 1:min(length(lev), length(typical_pressures))
        fprintf('  %2d  |       %3.0f\n', i, typical_pressures(i));
        if abs(typical_pressures(i) - 300) < 50  % Close to 300 hPa
            fprintf('      *** Closest to 300 hPa ***\n');
        end
    end
    
catch
    fprintf('Could not read pressure field. Using level index estimation.\n');
    fprintf('\nFor typical CAM hybrid coordinates:\n');
    fprintf('Level 7-9 are usually around 300 hPa\n');
    fprintf('Recommended: Try level 8 first\n');
end

% Sample wind data to check dimensions
fprintf('\nChecking wind data dimensions...\n');
try
    u_sample = ncread(hist_file, 'U', [1, 1, 1, 1], [1, 1, 1, 1]);
    fprintf('U wind data dimensions: time, lev, lat, lon\n');
    fprintf('Reading data at level 8 for testing...\n');
    
    u_test = squeeze(ncread(hist_file, 'U', [1, 8, 1, 1], [10, 1, 10, 10]));
    fprintf('Sample U wind values at level 8:\n');
    fprintf('  Min: %.2f m/s, Max: %.2f m/s, Mean: %.2f m/s\n', min(u_test(:)), max(u_test(:)), mean(u_test(:)));
    
    if all(isfinite(u_test(:))) && max(abs(u_test(:))) < 200
        fprintf('Level 8 data looks reasonable for 300 hPa winds.\n');
    else
        fprintf('Level 8 data may not be appropriate. Check other levels.\n');
    end
    
catch ME
    fprintf('Error reading wind data: %s\n', ME.message);
end

% Recommendations
fprintf('\n=== RECOMMENDATIONS ===\n');
fprintf('1. Based on typical CAM output, level 8 is likely close to 300 hPa\n');
fprintf('2. If you have access to the model documentation, verify the hybrid coordinate levels\n');
fprintf('3. You can also check levels 7 and 9 if level 8 doesn\'t look right\n');
fprintf('4. In the main script, set: lev_300hpa = 8; (or adjust as needed)\n');

fprintf('\nPress any key to continue...\n');
pause;