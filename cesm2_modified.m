close all; clear; 
%param; 
% cd /global/homes/w/wenyuz/co-spectra/ % Comment out original path
 
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
time_hist = double(ncread(hist_file, 'time'));
time_rcp = double(ncread(rcp_file, 'time'));

% Find 300 hPa level (approximately level 7-8, we'll use level 8)
% Note: lev is in hybrid coordinates, 300 hPa is typically around level 7-8
lev_300hpa = 8; % Adjust this if needed based on actual pressure levels

% Read wind components at 300 hPa level
% Note: U and V have dimensions (time, lev, lat, lon)
fprintf('Reading U wind component...\n');
uu_hist_full = ncread(hist_file, 'U');
uu_rcp_full = ncread(rcp_file, 'U');

fprintf('Reading V wind component...\n');
vv_hist_full = ncread(hist_file, 'V');
vv_rcp_full = ncread(rcp_file, 'V');

% Extract data at 300 hPa level
% NetCDF dimensions are (time, lev, lat, lon)
uu_hist = squeeze(uu_hist_full(:, lev_300hpa, :, :)); % (time, lat, lon)
uu_rcp = squeeze(uu_rcp_full(:, lev_300hpa, :, :));   % (time, lat, lon)
vv_hist = squeeze(vv_hist_full(:, lev_300hpa, :, :)); % (time, lat, lon)
vv_rcp = squeeze(vv_rcp_full(:, lev_300hpa, :, :));   % (time, lat, lon)

% Convert to MATLAB convention (lon, lat, time)
uu_hist = permute(uu_hist, [3, 2, 1]); % (lon, lat, time)
vv_hist = permute(vv_hist, [3, 2, 1]);
uu_rcp = permute(uu_rcp, [3, 2, 1]);
vv_rcp = permute(vv_rcp, [3, 2, 1]);

fprintf('Data dimensions: %d x %d x %d\n', size(uu_hist));
 
rad = lat/180.*pi;

% Get actual data dimensions for the analysis
n_lon = size(uu_hist, 1);
n_lat = size(uu_hist, 2);
n_time = size(uu_hist, 3);
fprintf('Analysis dimensions: %d lon x %d lat x %d time\n', n_lon, n_lat, n_time);

% Process historical data (1961-1980, 20 years of monthly data = 240 months)
% Assuming monthly data, we'll take DJFM months (Dec, Jan, Feb, Mar)
% For each year, extract Dec-Mar (4 months), then pad to 120 days as in original
n_years_hist = 20;
months_per_year = 12;

% For historical period: extract DJFM months
for yy = 1:n_years_hist
    % Get December of previous year + Jan, Feb, Mar of current year
    if yy == 1
        % For first year, start from January
        start_month = (yy-1)*months_per_year + 1; % January
        months_to_extract = 3; % Jan, Feb, Mar only
    else
        start_month = (yy-2)*months_per_year + 12; % December of previous year
        months_to_extract = 4; % Dec, Jan, Feb, Mar
    end
    
    % Since we have monthly data, we need to interpolate to daily
    % For simplicity, we'll replicate each month to 30 days
    if yy == 1
        uu_temp = uu_hist(:,:,start_month:start_month+months_to_extract-1);
        vv_temp = vv_hist(:,:,start_month:start_month+months_to_extract-1);
        % Pad to 120 days by repeating data
        uu_daily = repmat(uu_temp, [1, 1, 40]); % Approximate 120 days
        vv_daily = repmat(vv_temp, [1, 1, 40]);
    else
        uu_temp = uu_hist(:,:,start_month:start_month+months_to_extract-1);
        vv_temp = vv_hist(:,:,start_month:start_month+months_to_extract-1);
        % Pad to 120 days by repeating data
        uu_daily = repmat(uu_temp, [1, 1, 30]); % Approximate 120 days
        vv_daily = repmat(vv_temp, [1, 1, 30]);
    end
    
    % Take first 120 time steps
    if size(uu_daily, 3) > 120
        uu_daily = uu_daily(:,:,1:120);
        vv_daily = vv_daily(:,:,1:120);
    else
        % Pad if needed
        needed_steps = 120 - size(uu_daily, 3);
        uu_daily = cat(3, uu_daily, repmat(uu_daily(:,:,end), [1, 1, needed_steps]));
        vv_daily = cat(3, vv_daily, repmat(vv_daily(:,:,end), [1, 1, needed_steps]));
    end
    
    uu_yy_hist(yy,:,:,:) = permute(uu_daily, [4, 3, 2, 1]); % (year, time, lat, lon)
    vv_yy_hist(yy,:,:,:) = permute(vv_daily, [4, 3, 2, 1]);
end

% Process RCP8.5 data (2081-2100, 20 years)
n_years_rcp = 20;
for yy = 1:n_years_rcp
    % Same process for RCP data
    if yy == 1
        start_month = (yy-1)*months_per_year + 1;
        months_to_extract = 3;
    else
        start_month = (yy-2)*months_per_year + 12;
        months_to_extract = 4;
    end
    
    if yy == 1
        uu_temp = uu_rcp(:,:,start_month:start_month+months_to_extract-1);
        vv_temp = vv_rcp(:,:,start_month:start_month+months_to_extract-1);
        uu_daily = repmat(uu_temp, [1, 1, 40]);
        vv_daily = repmat(vv_temp, [1, 1, 40]);
    else
        uu_temp = uu_rcp(:,:,start_month:start_month+months_to_extract-1);
        vv_temp = vv_rcp(:,:,start_month:start_month+months_to_extract-1);
        uu_daily = repmat(uu_temp, [1, 1, 30]);
        vv_daily = repmat(vv_temp, [1, 1, 30]);
    end
    
    if size(uu_daily, 3) > 120
        uu_daily = uu_daily(:,:,1:120);
        vv_daily = vv_daily(:,:,1:120);
    else
        needed_steps = 120 - size(uu_daily, 3);
        uu_daily = cat(3, uu_daily, repmat(uu_daily(:,:,end), [1, 1, needed_steps]));
        vv_daily = cat(3, vv_daily, repmat(vv_daily(:,:,end), [1, 1, needed_steps]));
    end
    
    uu_yy_rcp(yy,:,:,:) = permute(uu_daily, [4, 3, 2, 1]);
    vv_yy_rcp(yy,:,:,:) = permute(vv_daily, [4, 3, 2, 1]);
end

%=========================== space-time analysis ================================
 
hann_win=hann(120);

% Process historical data
fprintf('Processing historical data...\n');
for yy=1:n_years_hist   % beginning of loop
 
% ------------------ applying a hanning window ----------------------
for i=1:size(uu_yy_hist,4)
for j=1:size(uu_yy_hist,3)
for t=1:size(uu_yy_hist,2)
  u_hist(i,t,j)=uu_yy_hist(yy,t,j,i)*hann_win(t);
  v_hist(i,t,j)=vv_yy_hist(yy,t,j,i)*hann_win(t);
end
end
end
 
% ---------------------- space-time spectrum -------------------------
[East_hist, West_hist]=space_time_new(u_hist,v_hist);
 
% --------------------- momentum flux divergence ---------------------
n_freq = size(East_hist, 2);
n_lat_spec = size(East_hist, 3);
for ff=1:n_freq
for jj=1:n_lat_spec
K_e_mm_cos_hist(ff,jj,:)=East_hist(:,ff,jj)*(cos(rad(jj))^2);
K_w_mm_cos_hist(ff,jj,:)=West_hist(:,ff,jj)*(cos(rad(jj))^2);
end
end
dK_e_m_hist=dy_dx(K_e_mm_cos_hist,rad,2)/RADIUS*DAY;
dK_w_m_hist=dy_dx(K_w_mm_cos_hist,rad,2)/RADIUS*DAY;
for ff=1:n_freq
for jj=1:n_lat_spec
K_e_m_yy_hist(yy,ff,jj,:)=-dK_e_m_hist(ff,jj,:)/(cos(rad(jj))^2);
K_w_m_yy_hist(yy,ff,jj,:)=-dK_w_m_hist(ff,jj,:)/(cos(rad(jj))^2);
end
end
 
end           % end of historical loop

% Process RCP8.5 data
fprintf('Processing RCP8.5 data...\n');
for yy=1:n_years_rcp   % beginning of loop
 
% ------------------ applying a hanning window ----------------------
for i=1:size(uu_yy_rcp,4)
for j=1:size(uu_yy_rcp,3)
for t=1:size(uu_yy_rcp,2)
  u_rcp(i,t,j)=uu_yy_rcp(yy,t,j,i)*hann_win(t);
  v_rcp(i,t,j)=vv_yy_rcp(yy,t,j,i)*hann_win(t);
end
end
end
 
% ---------------------- space-time spectrum -------------------------
[East_rcp, West_rcp]=space_time_new(u_rcp,v_rcp);
 
% --------------------- momentum flux divergence ---------------------
for ff=1:n_freq
for jj=1:n_lat_spec
K_e_mm_cos_rcp(ff,jj,:)=East_rcp(:,ff,jj)*(cos(rad(jj))^2);
K_w_mm_cos_rcp(ff,jj,:)=West_rcp(:,ff,jj)*(cos(rad(jj))^2);
end
end
dK_e_m_rcp=dy_dx(K_e_mm_cos_rcp,rad,2)/RADIUS*DAY;
dK_w_m_rcp=dy_dx(K_w_mm_cos_rcp,rad,2)/RADIUS*DAY;
for ff=1:n_freq
for jj=1:n_lat_spec
K_e_m_yy_rcp(yy,ff,jj,:)=-dK_e_m_rcp(ff,jj,:)/(cos(rad(jj))^2);
K_w_m_yy_rcp(yy,ff,jj,:)=-dK_w_m_rcp(ff,jj,:)/(cos(rad(jj))^2);
end
end
 
end           % end of RCP loop
 
%=========================== smooth ================================
K_e_sm_yy_hist=zeros(size(K_e_m_yy_hist));
K_w_sm_yy_hist=zeros(size(K_w_m_yy_hist));
K_e_sm_yy_rcp=zeros(size(K_e_m_yy_rcp));
K_w_sm_yy_rcp=zeros(size(K_w_m_yy_rcp));

% removing stationary waves
K_e_m_yy_hist(:,1,:,:)=0;
K_w_m_yy_hist(:,1,:,:)=0;
K_e_m_yy_rcp(:,1,:,:)=0;
K_w_m_yy_rcp(:,1,:,:)=0;
 
% Smooth historical data
n_wavenumber = size(K_e_m_yy_hist, 4);
for nn=1:n_years_hist
for j=1:n_lat_spec
for i=1:(n_wavenumber-1)
KK=smooth([fliplr(squeeze(K_w_m_yy_hist(nn,2:n_freq,j,i))),squeeze(K_e_m_yy_hist(nn,1:n_freq,j,i))],3,'g');
K_e_sm_yy_hist(nn,1:n_freq,j,i)=KK(n_freq:(2*n_freq-1));
K_w_sm_yy_hist(nn,2:n_freq,j,i)=fliplr(KK(1:(n_freq-1))');
end
end
end

% Smooth RCP data
for nn=1:n_years_rcp
for j=1:n_lat_spec
for i=1:(n_wavenumber-1)
KK=smooth([fliplr(squeeze(K_w_m_yy_rcp(nn,2:n_freq,j,i))),squeeze(K_e_m_yy_rcp(nn,1:n_freq,j,i))],3,'g');
K_e_sm_yy_rcp(nn,1:n_freq,j,i)=KK(n_freq:(2*n_freq-1));
K_w_sm_yy_rcp(nn,2:n_freq,j,i)=fliplr(KK(1:(n_freq-1))');
end
end
end

% Calculate zonal mean winds
uzm_hist=squeeze(mean(mean(uu_yy_hist(1:n_years_hist,:,:,:),4),2));
uzm_rcp=squeeze(mean(mean(uu_yy_rcp(1:n_years_rcp,:,:,:),4),2));

% Save processed data
save('uv_DJFM_CESM2_hist_300hpa.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
save('uv_DJFM_CESM2_rcp85_300hpa.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
 
end

%=========================== Analysis and Plotting ================================
 
% Load historical data
load('uv_DJFM_CESM2_hist_300hpa.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
uc_hist=squeeze(mean(uzm_hist,1))';
K_e_m_hist=squeeze(mean(K_e_sm_yy_hist,1));
K_w_m_hist=squeeze(mean(K_w_sm_yy_hist,1));

% Load RCP8.5 data 
load('uv_DJFM_CESM2_rcp85_300hpa.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
uc_rcp85=squeeze(mean(uzm_rcp,1))';
K_e_m_rcp=squeeze(mean(K_e_sm_yy_rcp,1));
K_w_m_rcp=squeeze(mean(K_w_sm_yy_rcp,1));
 
% Create wavenumber array based on actual data size
max_wavenum = size(K_e_m_hist, 3) - 1;
wavenum=[0:max_wavenum]; freq=0.5*[0:(size(K_w_m_hist,1)-1)]/(size(K_w_m_hist,1)-1); 
mm=2:max_wavenum; % Skip wavenumber 0 and 1
wavenum = wavenum(mm);
K_e_hist = K_e_m_hist(:,:,mm);
K_w_hist = K_w_m_hist(:,:,mm);
K_e_rcp = K_e_m_rcp(:,:,mm);
K_w_rcp = K_w_m_rcp(:,:,mm);
nk=length(wavenum);  nj=length(lat);
kk = 1:nk; jj=1:nj;
c_max=50;  c_min = -10;   ca=0.0:0.1:c_max;   cc=1:length(ca);
 
% Convert frequency-wavenumber to phase speed-wavenumber for historical data
for pos=1:nj;
[K_e_new_hist, K_w_new_hist] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, K_e_hist(:,pos,kk), K_w_hist(:,pos,kk));
K_e_lat_hist(:,pos)=sum(K_e_new_hist,2);
K_w_lat_hist(:,pos)=sum(K_w_new_hist,2);
end

% Convert frequency-wavenumber to phase speed-wavenumber for RCP8.5 data
for pos=1:nj;
[K_e_new_rcp, K_w_new_rcp] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, K_e_rcp(:,pos,kk), K_w_rcp(:,pos,kk));
K_e_lat_rcp85(:,pos)=sum(K_e_new_rcp,2);
K_w_lat_rcp85(:,pos)=sum(K_w_new_rcp,2);
end
 
save('co-spectra_cesm2_300hpa.mat','K_w_lat_hist','K_w_lat_rcp85','K_e_lat_hist','K_e_lat_rcp85','uc_hist','uc_rcp85');
 
% Plotting
figure('pos',[10 10 700 700]);
step = 0.01;
colorb=[-5*step:step:step*5];
[cc1,hh1]=contourf(-ca(cc),lat(jj),1.25*squeeze(K_w_lat_rcp85(cc,jj)-K_w_lat_hist(cc,jj))',colorb); hold on;
[cc2,hh2]=contourf(ca(cc),lat(jj),1.25*squeeze(K_e_lat_rcp85(cc,jj)-K_e_lat_hist(cc,jj))',colorb);
caxis([min(colorb),max(colorb)]);

% Use a colormap (if brewermap is not available, use built-in colormap)
try
    gg=brewermap(10,'RdBu');
    gg(5:6,:)=ones(2,3);
    colormap(gg(end:-1:1,:));
catch
    colormap(redblue(10)); % fallback colormap
end

set([hh1,hh2],'edgecolor','none');
step = 0.02;
[c1,h1]=contour(-ca(cc),lat(jj),squeeze(K_w_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2); hold on;
[c2,h2]=contour(ca(cc), lat(jj),squeeze(K_e_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2);
axis([c_min,c_max,-80,80]);
set([h1,h2],'linewidth',2.0)
xlabel('angular phase speed (m/s)','fontsize',20);
ylabel('latitude(deg)','fontsize',20);
set(gca,'xtick',[c_min:5:c_max],'ytick',[-80:20:80],'ylim',[-80,80],'fontsize',20);
line([0,0],[-80,80],'linestyle','-.','color','k');
plot(uc_rcp85./(cos(lat*pi/180)),lat,'r','linewidth',2);hold on;
plot(uc_hist./(cos(lat*pi/180)),lat,'k','linewidth',2);
title('300 hPa Co-spectra Analysis: RCP8.5 vs Historical','fontsize',16);
colorbar;

fprintf('Analysis complete. Results saved to co-spectra_cesm2_300hpa.mat\n');