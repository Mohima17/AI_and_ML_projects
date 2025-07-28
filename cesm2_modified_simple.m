close all; clear; 
%param; 
 
RADIUS=6.371e6;
DAY=86400;
 
do_compute=1; % Set to 1 to compute data
 
if do_compute>0 
%=========================== read data ================================
% Updated file paths and variable names for the provided NetCDF files
hist_file = 'b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc';
rcp_file = 'b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc';

% Read coordinate variables
lat = double(ncread(hist_file, 'lat'));
lon = double(ncread(hist_file, 'lon'));
lev = double(ncread(hist_file, 'lev'));

fprintf('Data grid: %d lon x %d lat x %d levels\n', length(lon), length(lat), length(lev));

% Find 300 hPa level (approximately level 7-8, we'll use level 8)
lev_300hpa = 8; % Adjust this if needed based on actual pressure levels

% Read wind components at 300 hPa level
fprintf('Reading wind data at level %d (300 hPa)...\n', lev_300hpa);
uu_hist_full = ncread(hist_file, 'U');
uu_rcp_full = ncread(rcp_file, 'U');
vv_hist_full = ncread(hist_file, 'V');
vv_rcp_full = ncread(rcp_file, 'V');

% Extract data at 300 hPa level
% NetCDF dimensions are (time, lev, lat, lon)
uu_hist = squeeze(uu_hist_full(:, lev_300hpa, :, :)); % (time, lat, lon)
uu_rcp = squeeze(uu_rcp_full(:, lev_300hpa, :, :));   
vv_hist = squeeze(vv_hist_full(:, lev_300hpa, :, :)); 
vv_rcp = squeeze(vv_rcp_full(:, lev_300hpa, :, :));   

fprintf('Wind data dimensions: %d time x %d lat x %d lon\n', size(uu_hist));

% Convert to MATLAB convention (lon, lat, time)
uu_hist = permute(uu_hist, [3, 2, 1]); % (lon, lat, time)
vv_hist = permute(vv_hist, [3, 2, 1]);
uu_rcp = permute(uu_rcp, [3, 2, 1]);
vv_rcp = permute(vv_rcp, [3, 2, 1]);

rad = lat/180.*pi;
n_lon = length(lon);
n_lat = length(lat);
n_time_hist = size(uu_hist, 3);
n_time_rcp = size(uu_rcp, 3);

fprintf('Final data structure: %d lon x %d lat x %d time (hist) x %d time (rcp)\n', ...
        n_lon, n_lat, n_time_hist, n_time_rcp);

%=========================== Extract seasonal data ================================
% For monthly data, extract DJF months (Dec, Jan, Feb) for each year
% Assuming data starts from January

% Process historical data (20 years, 240 months)
n_years_hist = 20;
djf_months_hist = [];
for yr = 1:n_years_hist
    start_idx = (yr-1)*12 + 1; % January of year yr
    if yr == 1
        % First year: only Jan, Feb (no previous December)
        djf_idx = [start_idx, start_idx+1]; % Jan, Feb
    else
        % Other years: Dec of previous year + Jan, Feb of current year
        djf_idx = [start_idx-1, start_idx, start_idx+1]; % Dec, Jan, Feb
    end
    djf_months_hist = [djf_months_hist, djf_idx];
end

% Process RCP data (20 years, 240 months)
n_years_rcp = 20;
djf_months_rcp = [];
for yr = 1:n_years_rcp
    start_idx = (yr-1)*12 + 1;
    if yr == 1
        djf_idx = [start_idx, start_idx+1];
    else
        djf_idx = [start_idx-1, start_idx, start_idx+1];
    end
    djf_months_rcp = [djf_months_rcp, djf_idx];
end

% Extract DJF data
uu_djf_hist = uu_hist(:, :, djf_months_hist);
vv_djf_hist = vv_hist(:, :, djf_months_hist);
uu_djf_rcp = uu_rcp(:, :, djf_months_rcp);
vv_djf_rcp = vv_rcp(:, :, djf_months_rcp);

fprintf('DJF data extracted: %d months (hist), %d months (rcp)\n', ...
        length(djf_months_hist), length(djf_months_rcp));

%=========================== Simple space-time analysis ================================
fprintf('Performing space-time analysis...\n');

% Apply window function
n_djf_hist = size(uu_djf_hist, 3);
n_djf_rcp = size(uu_djf_rcp, 3);

% Use available time steps (don't force to 120)
if n_djf_hist > 2
    hann_win_hist = hann(n_djf_hist);
    for t = 1:n_djf_hist
        uu_djf_hist(:, :, t) = uu_djf_hist(:, :, t) * hann_win_hist(t);
        vv_djf_hist(:, :, t) = vv_djf_hist(:, :, t) * hann_win_hist(t);
    end
end

if n_djf_rcp > 2
    hann_win_rcp = hann(n_djf_rcp);
    for t = 1:n_djf_rcp
        uu_djf_rcp(:, :, t) = uu_djf_rcp(:, :, t) * hann_win_rcp(t);
        vv_djf_rcp(:, :, t) = vv_djf_rcp(:, :, t) * hann_win_rcp(t);
    end
end

% Compute space-time spectra
fprintf('Computing space-time spectra for historical period...\n');
[East_hist, West_hist] = space_time_new(uu_djf_hist, vv_djf_hist);

fprintf('Computing space-time spectra for RCP period...\n');
[East_rcp, West_rcp] = space_time_new(uu_djf_rcp, vv_djf_rcp);

fprintf('Historical spectra dimensions: %d x %d x %d\n', size(East_hist));
fprintf('RCP spectra dimensions: %d x %d x %d\n', size(East_rcp));

%=========================== Momentum flux calculation ================================
fprintf('Computing momentum flux divergence...\n');

% Make sure both spectra have the same dimensions
min_freq = min(size(East_hist, 2), size(East_rcp, 2));
min_lat_spec = min(size(East_hist, 3), size(East_rcp, 3));

East_hist = East_hist(:, 1:min_freq, 1:min_lat_spec);
West_hist = West_hist(:, 1:min_freq, 1:min_lat_spec);
East_rcp = East_rcp(:, 1:min_freq, 1:min_lat_spec);
West_rcp = West_rcp(:, 1:min_freq, 1:min_lat_spec);

% Select latitude subset if needed (to match spectral output)
lat_spec = lat(1:min_lat_spec);
rad_spec = lat_spec/180.*pi;

% Compute momentum flux with cosine weighting
for ff = 1:min_freq
    for jj = 1:min_lat_spec
        K_e_cos_hist(ff, jj, :) = East_hist(:, ff, jj) * (cos(rad_spec(jj))^2);
        K_w_cos_hist(ff, jj, :) = West_hist(:, ff, jj) * (cos(rad_spec(jj))^2);
        K_e_cos_rcp(ff, jj, :) = East_rcp(:, ff, jj) * (cos(rad_spec(jj))^2);
        K_w_cos_rcp(ff, jj, :) = West_rcp(:, ff, jj) * (cos(rad_spec(jj))^2);
    end
end

% Compute meridional derivative
dK_e_hist = dy_dx(K_e_cos_hist, rad_spec, 2) / RADIUS * DAY;
dK_w_hist = dy_dx(K_w_cos_hist, rad_spec, 2) / RADIUS * DAY;
dK_e_rcp = dy_dx(K_e_cos_rcp, rad_spec, 2) / RADIUS * DAY;
dK_w_rcp = dy_dx(K_w_cos_rcp, rad_spec, 2) / RADIUS * DAY;

% Final momentum flux divergence
for ff = 1:min_freq
    for jj = 1:min_lat_spec
        K_e_flux_hist(ff, jj, :) = -dK_e_hist(ff, jj, :) / (cos(rad_spec(jj))^2);
        K_w_flux_hist(ff, jj, :) = -dK_w_hist(ff, jj, :) / (cos(rad_spec(jj))^2);
        K_e_flux_rcp(ff, jj, :) = -dK_e_rcp(ff, jj, :) / (cos(rad_spec(jj))^2);
        K_w_flux_rcp(ff, jj, :) = -dK_w_rcp(ff, jj, :) / (cos(rad_spec(jj))^2);
    end
end

% Smooth the results (simplified)
K_e_smooth_hist = K_e_flux_hist;
K_w_smooth_hist = K_w_flux_hist;
K_e_smooth_rcp = K_e_flux_rcp;
K_w_smooth_rcp = K_w_flux_rcp;

% Remove stationary waves (wavenumber 0)
K_e_smooth_hist(1, :, :) = 0;
K_w_smooth_hist(1, :, :) = 0;
K_e_smooth_rcp(1, :, :) = 0;
K_w_smooth_rcp(1, :, :) = 0;

% Calculate zonal mean winds
uzm_hist = squeeze(mean(uu_djf_hist, 1)); % Average over longitude
uzm_rcp = squeeze(mean(uu_djf_rcp, 1));
uc_hist = squeeze(mean(uzm_hist, 2))'; % Average over time
uc_rcp = squeeze(mean(uzm_rcp, 2))';

% Save results
save('co-spectra_cesm2_300hpa_simple.mat', 'K_e_smooth_hist', 'K_w_smooth_hist', ...
     'K_e_smooth_rcp', 'K_w_smooth_rcp', 'uc_hist', 'uc_rcp', 'lat_spec');

fprintf('Data processing complete. Results saved.\n');

end

%=========================== Analysis and Plotting ================================

% Load data
load('co-spectra_cesm2_300hpa_simple.mat');

% Create wavenumber and frequency arrays
max_wavenum = size(K_e_smooth_hist, 3) - 1;
wavenum = 2:max_wavenum; % Skip wavenumbers 0 and 1
freq = 0.5 * [0:(size(K_e_smooth_hist, 1)-1)] / (size(K_e_smooth_hist, 1)-1);

% Select subset for analysis
K_e_hist_sub = K_e_smooth_hist(:, :, wavenum);
K_w_hist_sub = K_w_smooth_hist(:, :, wavenum);
K_e_rcp_sub = K_e_smooth_rcp(:, :, wavenum);
K_w_rcp_sub = K_w_smooth_rcp(:, :, wavenum);

% Convert to phase speed representation (simplified)
c_max = 50; c_min = -10; ca = 0.0:0.1:c_max; cc = 1:length(ca);
nj = length(lat_spec); jj = 1:nj;

% Sum over wavenumbers for each latitude (simplified approach)
K_e_lat_hist = squeeze(sum(K_e_hist_sub, 3));
K_w_lat_hist = squeeze(sum(K_w_hist_sub, 3));
K_e_lat_rcp = squeeze(sum(K_e_rcp_sub, 3));
K_w_lat_rcp = squeeze(sum(K_w_rcp_sub, 3));

% Basic plotting
figure('Position', [10 10 700 700]);
step = 0.01;
colorb = [-5*step:step:step*5];

% Create simple difference plot
lat_plot = lat_spec(jj);
freq_plot = freq;

% Plot differences (simplified)
subplot(2,1,1);
contourf(freq_plot, lat_plot, (K_e_lat_rcp - K_e_lat_hist)', colorb);
colorbar;
title('Eastward Wave Changes (RCP8.5 - Historical) at 300 hPa');
xlabel('Frequency');
ylabel('Latitude (deg)');

subplot(2,1,2);
contourf(freq_plot, lat_plot, (K_w_lat_rcp - K_w_lat_hist)', colorb);
colorbar;
title('Westward Wave Changes (RCP8.5 - Historical) at 300 hPa');
xlabel('Frequency');
ylabel('Latitude (deg)');

% Plot zonal winds
figure;
plot(uc_hist, lat_spec, 'k-', 'LineWidth', 2); hold on;
plot(uc_rcp, lat_spec, 'r-', 'LineWidth', 2);
xlabel('Zonal Wind (m/s)');
ylabel('Latitude (deg)');
legend('Historical', 'RCP8.5', 'Location', 'best');
title('300 hPa Zonal Wind Comparison');
grid on;

fprintf('Analysis complete. Figures generated.\n');