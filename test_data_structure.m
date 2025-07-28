% Quick test script to verify data structure and find correct 300 hPa level
clear; close all;

% File names
hist_file = 'b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc';

if ~exist(hist_file, 'file')
    fprintf('Error: File %s not found in current directory\n', hist_file);
    fprintf('Current directory: %s\n', pwd);
    fprintf('Please ensure the NetCDF files are in the correct location.\n');
    return;
end

fprintf('Testing data structure...\n');

% Read basic info
try
    info = ncinfo(hist_file);
    fprintf('NetCDF file structure:\n');
    for i = 1:length(info.Variables)
        var = info.Variables(i);
        if length(var.Size) > 1
            fprintf('  %s: ', var.Name);
            fprintf('%d ', var.Size);
            fprintf('\n');
        end
    end
catch ME
    fprintf('Error reading file info: %s\n', ME.message);
    return;
end

% Read coordinate variables
try
    lat = double(ncread(hist_file, 'lat'));
    lon = double(ncread(hist_file, 'lon'));
    lev = double(ncread(hist_file, 'lev'));
    time = double(ncread(hist_file, 'time'));
    
    fprintf('\nCoordinate dimensions:\n');
    fprintf('  Longitude: %d points (%.1f to %.1f)\n', length(lon), min(lon), max(lon));
    fprintf('  Latitude: %d points (%.1f to %.1f)\n', length(lat), min(lat), max(lat));
    fprintf('  Levels: %d levels\n', length(lev));
    fprintf('  Time: %d time steps\n', length(time));
    
catch ME
    fprintf('Error reading coordinates: %s\n', ME.message);
    return;
end

% Test reading U wind data
try
    fprintf('\nTesting U wind data reading...\n');
    u_info = ncinfo(hist_file, 'U');
    fprintf('U wind dimensions: ');
    fprintf('%s ', u_info.Dimensions.Name);
    fprintf('\n');
    fprintf('U wind size: ');
    fprintf('%d ', u_info.Size);
    fprintf('\n');
    
    % Try reading a small sample
    u_sample = ncread(hist_file, 'U', [1, 1, 1, 1], [2, 2, 2, 2]);
    fprintf('Sample U data size: %d x %d x %d x %d\n', size(u_sample));
    fprintf('Sample U values: min=%.2f, max=%.2f\n', min(u_sample(:)), max(u_sample(:)));
    
catch ME
    fprintf('Error reading U wind: %s\n', ME.message);
    return;
end

% Test different level extractions
fprintf('\nTesting level extraction for 300 hPa...\n');
for test_lev = 7:10
    try
        if test_lev <= length(lev)
            u_test = squeeze(ncread(hist_file, 'U', [1, test_lev, 1, 1], [10, 1, 10, 10]));
            fprintf('Level %d: size %d x %d, range %.1f to %.1f m/s\n', ...
                    test_lev, size(u_test, 1), size(u_test, 2), min(u_test(:)), max(u_test(:)));
            
            % Check if values are reasonable for upper atmosphere
            if max(abs(u_test(:))) > 10 && max(abs(u_test(:))) < 200
                fprintf('  --> Level %d looks good for 300 hPa\n', test_lev);
            end
        end
    catch
        fprintf('Level %d: Error reading data\n', test_lev);
    end
end

% Recommend level
fprintf('\n=== RECOMMENDATIONS ===\n');
fprintf('1. Use level 8 or 9 for 300 hPa (typical values 20-100 m/s)\n');
fprintf('2. Data dimensions are (time, lev, lat, lon)\n');
fprintf('3. After extraction at one level: (time, lat, lon)\n');
fprintf('4. For analysis, permute to (lon, lat, time)\n');

fprintf('\nTest complete. You can now run the simplified analysis script.\n');